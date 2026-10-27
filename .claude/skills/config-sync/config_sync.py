# Config Sync - add, retire or rename one key across every INI config variant
#
# Edits the INI files as text, line by line.  It does NOT use
# configparser.write(), which would drop every comment, decorative rule and
# blank line.  Each file keeps its own line endings (CRLF or LF), its
# encoding and whether or not it ends with a newline.  Old values are never
# deleted: a retired key is commented out ("#key = value").
#
# Default is a dry run that prints the diff.  Add "-w True" to write.
#
# Contributed by CLB Engineering for Andy Carter, PE -- 2026.10.27
# Standard library only (Python 3.8 or later)
# ************************************************************

# ************************************************************
import os
import re
import sys
import glob
import difflib

import argparse
import time
import datetime
# ************************************************************


# ----------------
def fn_str_to_bool(value):
    if isinstance(value, bool):
        return value
    if value.lower() in {'true', 't', '1'}:
        return True
    elif value.lower() in {'false', 'f', '0'}:
        return False
    else:
        raise argparse.ArgumentTypeError(f"Boolean value expected. Got '{value}'.")
# ----------------


# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
def fn_read_ini(str_path):
    # Returns a dict describing the file so it can be written back byte-for-byte
    with open(str_path, 'rb') as f:
        bytes_data = f.read()
    return fn_parse_ini_bytes(str_path, bytes_data)
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
def fn_parse_ini_bytes(str_path, bytes_data):
    if bytes_data.startswith(b'\xef\xbb\xbf'):
        str_encoding = 'utf-8-sig'
    else:
        try:
            bytes_data.decode('utf-8')
            str_encoding = 'utf-8'
        except UnicodeDecodeError:
            str_encoding = 'latin-1'

    str_text = bytes_data.decode(str_encoding)

    int_crlf = str_text.count('\r\n')
    int_lf = str_text.count('\n') - int_crlf

    if int_crlf and int_lf:
        str_newline = 'MIXED'
    elif int_crlf:
        str_newline = '\r\n'
    else:
        str_newline = '\n'

    list_lines = []
    b_final_newline = False
    if str_newline != 'MIXED':
        list_lines = str_text.split(str_newline)
        if len(list_lines) > 1 and list_lines[-1] == '':
            b_final_newline = True
            list_lines = list_lines[:-1]

    return {'path': str_path,
            'bytes': bytes_data,
            'encoding': str_encoding,
            'newline': str_newline,
            'final_newline': b_final_newline,
            'lines': list_lines}
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
def fn_build_bytes(dict_ini, list_new_lines):
    str_text = dict_ini['newline'].join(list_new_lines)
    if dict_ini['final_newline']:
        str_text += dict_ini['newline']
    return str_text.encode(dict_ini['encoding'])
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


# ----------------
def fn_describe(dict_ini):
    if dict_ini['newline'] == '\r\n':
        str_nl = 'CRLF'
    elif dict_ini['newline'] == '\n':
        str_nl = 'LF'
    else:
        str_nl = 'MIXED'
    str_final = 'final newline' if dict_ini['final_newline'] else 'no final newline'
    return f"{str_nl}, {dict_ini['encoding']}, {str_final}"
# ----------------


# ......................
def fn_section_spans(list_lines):
    # dict of section name -> (header index, end index exclusive)
    re_header = re.compile(r'^\s*\[([^\]]+)\]\s*$')
    list_headers = []
    for int_idx, str_line in enumerate(list_lines):
        m = re_header.match(str_line)
        if m:
            list_headers.append((m.group(1), int_idx))

    dict_spans = {}
    for int_i, (str_name, int_start) in enumerate(list_headers):
        if int_i + 1 < len(list_headers):
            int_end = list_headers[int_i + 1][1]
        else:
            int_end = len(list_lines)
        if str_name not in dict_spans:
            dict_spans[str_name] = (int_start, int_end)
    return dict_spans
# ......................


# ......................
def fn_active_key(str_line):
    # key name if the line is an active "key = value" (or "key: value") line
    if not str_line.strip() or str_line[:1].isspace():
        return None
    if str_line.lstrip()[:1] in ('#', ';', '['):
        return None
    m = re.match(r'^([^=:]+?)\s*[=:]', str_line)
    return m.group(1).strip() if m else None
# ......................


# ......................
def fn_commented_key(str_line):
    m = re.match(r'^\s*[#;]+\s*([A-Za-z0-9_.\-]+)\s*[=:]', str_line)
    return m.group(1) if m else None
# ......................


# ......................
def fn_find_key(list_lines, tup_span, str_key):
    # index of the active line holding str_key inside the section, or None
    # (configparser keys are not case sensitive, so neither is this)
    for int_idx in range(tup_span[0] + 1, tup_span[1]):
        str_found = fn_active_key(list_lines[int_idx])
        if str_found is not None and str_found.lower() == str_key.lower():
            return int_idx
    return None
# ......................


# ......................
def fn_find_commented(list_lines, tup_span, str_key):
    list_idx = []
    for int_idx in range(tup_span[0] + 1, tup_span[1]):
        str_found = fn_commented_key(list_lines[int_idx])
        if str_found is not None and str_found.lower() == str_key.lower():
            list_idx.append(int_idx)
    return list_idx
# ......................


# ......................
def fn_end_of_value(list_lines, int_idx, int_end):
    # skip indented continuation lines of a multi-line value
    while (int_idx + 1 < int_end and list_lines[int_idx + 1].strip()
           and list_lines[int_idx + 1][:1].isspace()):
        int_idx += 1
    return int_idx
# ......................


# ......................
def fn_comment_lines(str_comment):
    if not str_comment:
        return []
    return ['# ' + s.strip() for s in str_comment.split('\\n')]
# ......................


# ----------------------------------------
def fn_change_add(list_lines, str_section, str_key, str_value, str_comment, str_after):
    dict_spans = fn_section_spans(list_lines)
    if str_section not in dict_spans:
        return None, f"SKIPPED - section [{str_section}] not in this file"
    tup_span = dict_spans[str_section]

    int_existing = fn_find_key(list_lines, tup_span, str_key)
    if int_existing is not None:
        return None, f"SKIPPED - already present (line {int_existing + 1}: {list_lines[int_existing].strip()})"

    if str_after:
        int_anchor = fn_find_key(list_lines, tup_span, str_after)
        if int_anchor is None:
            return None, f"SKIPPED - anchor key '{str_after}' not active in [{str_section}]"
    else:
        # after the last active key in the section, else right after the header
        int_anchor = tup_span[0]
        for int_idx in range(tup_span[0] + 1, tup_span[1]):
            if fn_active_key(list_lines[int_idx]) is not None:
                int_anchor = int_idx
    int_anchor = fn_end_of_value(list_lines, int_anchor, tup_span[1])

    list_block = [''] + fn_comment_lines(str_comment) + [f"{str_key} = {str_value}"]
    list_new = list_lines[:int_anchor + 1] + list_block + list_lines[int_anchor + 1:]

    str_note = ''
    if fn_find_commented(list_lines, tup_span, str_key):
        str_note = ' (a commented-out copy of this key already exists in the section)'
    return list_new, f"ADDED at line {int_anchor + 1 + len(list_block)}{str_note}"
# ----------------------------------------


# ----------------------------------------
def fn_change_retire(list_lines, str_section, str_key, str_comment, str_new_key):
    # str_new_key = None for retire; a name for rename (old line kept as comment)
    dict_spans = fn_section_spans(list_lines)
    if str_section not in dict_spans:
        return None, f"SKIPPED - section [{str_section}] not in this file"
    tup_span = dict_spans[str_section]

    int_idx = fn_find_key(list_lines, tup_span, str_key)
    if int_idx is None:
        if fn_find_commented(list_lines, tup_span, str_key):
            return None, "SKIPPED - key is already commented out"
        return None, "SKIPPED - key not in this section"

    if str_new_key and fn_find_key(list_lines, tup_span, str_new_key) is not None:
        return None, f"SKIPPED - new key '{str_new_key}' already present"

    int_last = fn_end_of_value(list_lines, int_idx, tup_span[1])
    list_old = list_lines[int_idx:int_last + 1]

    str_today = datetime.date.today().strftime('%Y.%m.%d')
    list_block = []
    if str_comment:
        str_word = 'Renamed' if str_new_key else 'Retired'
        list_block.append(f"# {str_word} - {str_today} -- {str_comment}")
    list_block += ['#' + s for s in list_old]

    if str_new_key:
        str_first = list_old[0]
        int_delim = re.search(r'[=:]', str_first).start()
        str_value = str_first[int_delim + 1:].strip()
        list_rest = list_old[1:]
        list_block.append(f"{str_new_key} = {str_value}")
        list_block += list_rest
        str_msg = f"RENAMED at line {int_idx + 1}: {str_key} -> {str_new_key} (old line kept as comment)"
    else:
        str_msg = f"RETIRED at line {int_idx + 1} (commented out, value kept)"

    list_new = list_lines[:int_idx] + list_block + list_lines[int_last + 1:]
    return list_new, str_msg
# ----------------------------------------


# ----------------------------------------
def fn_change_add_section(list_lines, str_section, list_block, str_after):
    dict_spans = fn_section_spans(list_lines)
    if str_section in dict_spans:
        return None, f"SKIPPED - section [{str_section}] already present"

    if str_after:
        if str_after not in dict_spans:
            return None, f"SKIPPED - anchor section [{str_after}] not in this file"
        int_end = dict_spans[str_after][1]
    else:
        int_end = len(list_lines)

    # comment / rule lines directly above the next header belong to the next section
    int_lead = int_end
    if int_end < len(list_lines):
        while int_lead - 1 >= 0 and list_lines[int_lead - 1].strip().startswith(('#', ';')):
            int_lead -= 1

    int_body = int_lead
    while int_body - 1 >= 0 and not list_lines[int_body - 1].strip():
        int_body -= 1

    # blank lines between sections: copy what this file already uses
    int_sep = int_lead - int_body
    if int_end >= len(list_lines):
        int_sep = fn_typical_gap(list_lines)
    int_sep = max(int_sep, 1)

    list_insert = [''] * int_sep + list_block
    if int_end < len(list_lines) and int_lead == int_body:
        list_insert = list_insert + [''] * int_sep

    list_new = list_lines[:int_body] + list_insert + list_lines[int_body:]
    return list_new, f"SECTION ADDED at line {int_body + int_sep + 1}"
# ----------------------------------------


# ......................
def fn_typical_gap(list_lines):
    # blank lines before the last section's lead-in, used when appending at the end
    dict_spans = fn_section_spans(list_lines)
    if len(dict_spans) < 2:
        return 1
    int_header = max(tup[0] for tup in dict_spans.values())
    int_lead = int_header
    while int_lead - 1 >= 0 and list_lines[int_lead - 1].strip().startswith(('#', ';')):
        int_lead -= 1
    int_gap = 0
    while int_lead - 1 - int_gap >= 0 and not list_lines[int_lead - 1 - int_gap].strip():
        int_gap += 1
    return int_gap
# ......................


# ++++++++++++++++++++++++++++
def fn_verify(dict_before, dict_after):
    # confirm only intended lines changed: insertions, and "#" added in front
    list_problems = []
    if dict_before['newline'] != dict_after['newline']:
        list_problems.append('line endings changed')
    if dict_before['final_newline'] != dict_after['final_newline']:
        list_problems.append('final newline changed')
    if dict_before['encoding'] != dict_after['encoding']:
        list_problems.append('encoding changed')

    sm = difflib.SequenceMatcher(None, dict_before['lines'], dict_after['lines'], autojunk=False)
    for str_tag, i1, i2, j1, j2 in sm.get_opcodes():
        if str_tag in ('equal', 'insert'):
            continue
        if str_tag == 'delete':
            list_problems.append(f"lines {i1 + 1}-{i2} deleted")
            continue
        # replace: every old line must survive as "#" + old line, in order
        list_old = dict_before['lines'][i1:i2]
        list_new = dict_after['lines'][j1:j2]
        int_pos = 0
        for str_line in list_new:
            if int_pos < len(list_old) and str_line == '#' + list_old[int_pos]:
                int_pos += 1
        if int_pos != len(list_old):
            list_problems.append(f"lines {i1 + 1}-{i2} altered")
    return list_problems
# ++++++++++++++++++++++++++++


# ----------------
def fn_show(dict_ini, str_section, str_key):
    list_lines = dict_ini['lines']
    dict_spans = fn_section_spans(list_lines)
    if str_section not in dict_spans:
        return f"[{str_section}] not present"
    tup_span = dict_spans[str_section]
    if not str_key:
        return f"[{str_section}] at line {tup_span[0] + 1}"
    int_idx = fn_find_key(list_lines, tup_span, str_key)
    str_out = f"[{str_section}] "
    if int_idx is None:
        str_out += f"{str_key}: not active"
    else:
        str_out += f"line {int_idx + 1}: {list_lines[int_idx].strip()}"
    list_commented = fn_find_commented(list_lines, tup_span, str_key)
    if list_commented:
        str_out += f"  (commented copies at lines {', '.join(str(i + 1) for i in list_commented)})"
    return str_out
# ----------------


# +++++++++++++++++++++++++++++
def fn_config_sync(str_action,
                   list_file_patterns,
                   str_section,
                   str_key,
                   str_value,
                   str_new_key,
                   str_comment,
                   str_after,
                   str_block_file,
                   b_write):

    print(" ")
    print("+=================================================================+")
    print("|          CONFIG SYNC - ONE KEY ACROSS EVERY INI VARIANT         |")
    print("|   Contributed by CLB Engineering for Andy Carter, PE (2026)     |")
    print("+-----------------------------------------------------------------+")
    print("  ---(t) ACTION: " + str_action)
    print("  ---(s) SECTION: [" + str_section + "]")
    if str_key:
        print("  ---(k) KEY: " + str_key)
    print("  ---[w] WRITE CHANGES: " + str(b_write))
    print("===================================================================")

    list_files = []
    for str_pattern in list_file_patterns:
        list_match = sorted(glob.glob(str_pattern))
        if not list_match:
            print(f"  !! no file matches: {str_pattern}")
        for str_path in list_match:
            if str_path not in list_files:
                list_files.append(str_path)

    list_block = []
    if str_action == 'section':
        if str_block_file == '-':
            list_block = sys.stdin.read().splitlines()
        else:
            with open(str_block_file, 'r', encoding='utf-8') as f:
                list_block = f.read().splitlines()
        while list_block and not list_block[-1].strip():
            list_block.pop()

    int_changed = 0
    int_failed = 0
    for str_path in list_files:
        dict_ini = fn_read_ini(str_path)
        str_name = os.path.basename(str_path)

        if dict_ini['newline'] == 'MIXED':
            print(f"  -- {str_name}: SKIPPED - mixed CRLF/LF line endings, edit by hand")
            continue

        if str_action == 'show':
            print(f"  -- {str_name} [{fn_describe(dict_ini)}]: "
                  + fn_show(dict_ini, str_section, str_key))
            continue

        if str_action == 'add':
            list_new, str_msg = fn_change_add(dict_ini['lines'], str_section, str_key,
                                              str_value, str_comment, str_after)
        elif str_action == 'retire':
            list_new, str_msg = fn_change_retire(dict_ini['lines'], str_section, str_key,
                                                 str_comment, None)
        elif str_action == 'rename':
            list_new, str_msg = fn_change_retire(dict_ini['lines'], str_section, str_key,
                                                 str_comment, str_new_key)
        else:
            list_new, str_msg = fn_change_add_section(dict_ini['lines'], str_section,
                                                      list_block, str_after)

        print(f"  -- {str_name} [{fn_describe(dict_ini)}]: {str_msg}")
        if list_new is None:
            continue

        for str_diff in difflib.unified_diff(dict_ini['lines'], list_new,
                                             fromfile=str_path, tofile=str_path,
                                             lineterm='', n=2):
            print("     " + str_diff)

        bytes_new = fn_build_bytes(dict_ini, list_new)
        dict_check = fn_parse_ini_bytes(str_path, bytes_new)
        list_problems = fn_verify(dict_ini, dict_check)
        if list_problems:
            print("     !! NOT WRITTEN - check failed: " + '; '.join(list_problems))
            int_failed += 1
            continue

        if b_write:
            with open(str_path, 'wb') as f:
                f.write(bytes_new)
            dict_after = fn_read_ini(str_path)
            list_problems = fn_verify(dict_ini, dict_after)
            if list_problems:
                print("     !! WRITTEN BUT CHECK FAILED: " + '; '.join(list_problems))
                int_failed += 1
            else:
                print("     written; line endings, encoding and final newline unchanged")
        int_changed += 1

    print("+-----------------------------------------------------------------+")
    if str_action != 'show':
        str_mode = 'written' if b_write else 'would change (dry run, add -w True to write)'
        print(f"  Files {str_mode}: {int_changed} of {len(list_files)}")
        if int_failed:
            print(f"  Files failing the byte check: {int_failed}")
# +++++++++++++++++++++++++++++


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == '__main__':

    flt_start_run = time.time()

    parser = argparse.ArgumentParser(description='========= CONFIG SYNC - ONE KEY ACROSS EVERY INI VARIANT =========')

    parser.add_argument('-t',
                        dest = "str_action",
                        help=r'REQUIRED: add, retire, rename, section (add a whole section) or show',
                        required=True,
                        choices=['add', 'retire', 'rename', 'section', 'show'],
                        metavar='ACTION',
                        type=str)

    parser.add_argument('-f',
                        dest = "list_file_patterns",
                        help=r'REQUIRED: INI files or glob patterns Example: "src/config_global*.ini"',
                        required=True,
                        nargs='+',
                        metavar='FILE',
                        type=str)

    parser.add_argument('-s',
                        dest = "str_section",
                        help=r'REQUIRED: Section name without brackets Example: 08_filter_netcdf',
                        required=True,
                        metavar='SECTION',
                        type=str)

    parser.add_argument('-k',
                        dest = "str_key",
                        help=r'OPTIONAL: Key to add, retire, rename or show Default: None',
                        required=False,
                        default=None,
                        metavar='KEY',
                        type=str)

    parser.add_argument('-v',
                        dest = "str_value",
                        help=r'OPTIONAL: Value for "add" (same value in every file) Default: None',
                        required=False,
                        default=None,
                        metavar='VALUE',
                        type=str)

    parser.add_argument('-n',
                        dest = "str_new_key",
                        help=r'OPTIONAL: New key name for "rename" Default: None',
                        required=False,
                        default=None,
                        metavar='KEY',
                        type=str)

    parser.add_argument('-m',
                        dest = "str_comment",
                        help=r'OPTIONAL: Comment written above the change; "\n" splits lines Default: None',
                        required=False,
                        default=None,
                        metavar='TEXT',
                        type=str)

    parser.add_argument('-a',
                        dest = "str_after",
                        help=r'OPTIONAL: For "add", the key to insert after; for "section", the section to insert after Default: end of section / end of file',
                        required=False,
                        default=None,
                        metavar='NAME',
                        type=str)

    parser.add_argument('-b',
                        dest = "str_block_file",
                        help=r'OPTIONAL: For "section", a text file holding the whole new section (rules, header, comments, keys), or - to read it from stdin Default: None',
                        required=False,
                        default=None,
                        metavar='FILE',
                        type=str)

    parser.add_argument('-w',
                        dest = "b_write",
                        help=r'OPTIONAL: Write the changes (otherwise dry run, diff only) Default: False',
                        required=False,
                        default=False,
                        metavar='T/F',
                        type=fn_str_to_bool)

    args = vars(parser.parse_args())

    if args['str_action'] in ('add', 'retire', 'rename') and not args['str_key']:
        parser.error('-k KEY is required for add, retire and rename')
    if args['str_action'] == 'add' and args['str_value'] is None:
        parser.error('-v VALUE is required for add')
    if args['str_action'] == 'rename' and not args['str_new_key']:
        parser.error('-n NEW_KEY is required for rename')
    if args['str_action'] == 'section' and not args['str_block_file']:
        parser.error('-b BLOCK_FILE is required for section')

    fn_config_sync(args['str_action'],
                   args['list_file_patterns'],
                   args['str_section'],
                   args['str_key'],
                   args['str_value'],
                   args['str_new_key'],
                   args['str_comment'],
                   args['str_after'],
                   args['str_block_file'],
                   args['b_write'])

    flt_end_run = time.time()
    flt_time_pass = (flt_end_run - flt_start_run) // 1
    time_pass = datetime.timedelta(seconds=flt_time_pass)

    print('Compute Time: ' + str(time_pass))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
