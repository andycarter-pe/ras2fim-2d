<!-- Contributed by CLB Engineering for Andy Carter, PE, 2026.10.27. Supporting file for /run-demo. -->

# /run-demo reference (ras2fim-2d)

## Download the sample data without the AWS CLI
The bucket `rasfim-2d-sample` is public and holds fewer than 1,000 files, so one listing call is
enough. Set `PREFIX` to `sample_input/` or `sample_output/` and `DEST` to the target folder.

```
python - <<'EOF'
import os, re, urllib.parse, urllib.request
PREFIX, DEST = "sample_input/", "ras2fim2d_sample"
str_url = "https://rasfim-2d-sample.s3.amazonaws.com/"
str_xml = urllib.request.urlopen(str_url + "?list-type=2&prefix=" + PREFIX).read().decode()
for str_key in re.findall(r"<Key>(.*?)</Key>", str_xml):
    if str_key.endswith("/"):
        continue
    str_out = os.path.join(DEST, *str_key.split("/"))
    os.makedirs(os.path.dirname(str_out), exist_ok=True)
    urllib.request.urlretrieve(str_url + urllib.parse.quote(str_key), str_out)
    print(str_key)
EOF
```

One key in `sample_output/02b_prep_for_ras/` has a Windows path (`E:\...dss`) as its name. It is
a stray copy of a DSS file and is not needed by any step.

## Read NetCDF attributes
Run with the image's Python (`docker run --rm -v <run_folder>:/model_output <image> python - <<'EOF'`)
or the local env. Paths are inside the run folder.

```
import glob, netCDF4
for str_nc in sorted(glob.glob('/model_output/06_simple_rasters/0[23]_wsel_nc*/*.nc')):
    with netCDF4.Dataset(str_nc) as ds:
        print('==', str_nc)
        print('   dims:', {k: len(v) for k, v in ds.dimensions.items()})
        print('   vars:', list(ds.variables))
        for str_att in ds.ncattrs():
            print('   ', str_att, '=', str(ds.getncattr(str_att))[:200])
        if 'flow' in ds.variables:
            arr_flow = ds.variables['flow'][:]
            print('   flows:', len(arr_flow), 'from', arr_flow.min(), 'to', arr_flow.max(),
                  ds.variables['flow'].getncattr('units') if 'units' in ds.variables['flow'].ncattrs() else '')
```

## Expected sample results (from the S3 sample_output and the repo's sample_output)
- HEC-RAS run folder: `1919912_wb-2427466_wb-2427467_14-hr_50-cfs_to_11609-cfs` (the 50 cfs low
  flow comes from `flt_min_q = 50` in `config_global_windows.ini`; `config_global.ini` has 100).
- Reaches: `wb-2427466` and `wb-2427467`.
- Final products: `06_simple_rasters/03_wsel_nc_filtered/filter_wsel_wb-2427466.nc` and
  `filter_wsel_wb-2427467.nc`.
