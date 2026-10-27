"""Write a self-contained pixel review of one healed-terrain window."""

import html
import json
from pathlib import Path

import numpy as np


def write_review(model_path, output, qa):
    """Compare aligned original/healed/difference pixels in the first raised block.

    Reads at most 256 square cells per preview. This is a locator for manual GIS
    review, not a full-extent map or engineering acceptance. Inputs are read-only.
    """
    from terrio import rio
    from terrio.rio.windows import Window

    output = Path(output)
    with rio.open(output / "validity.tif") as mask:
        window = Window(0, 0, min(256, mask.width), min(256, mask.height))
        found = False
        for row in range(0, mask.height, 256):
            for col in range(0, mask.width, 256):
                candidate = Window(
                    col, row, min(256, mask.width - col), min(256, mask.height - row)
                )
                if mask.read(1, window=candidate).any():
                    window, found = candidate, True
                    break
            if found:
                break
    arrays = []
    for path in (model_path, output / "healed.tif", output / "difference.tif"):
        with rio.open(path) as src:
            values = src.read(1, window=window).astype(float)
            valid = np.isfinite(values) & (values != src.nodata)
            arrays.append(np.where(valid, values, np.nan))
    elevations = np.concatenate([a[np.isfinite(a)] for a in arrays[:2]])
    low, high = (
        (float(elevations.min()), float(elevations.max()))
        if elevations.size
        else (0.0, 1.0)
    )
    change = arrays[2][np.isfinite(arrays[2])]
    maximum = float(change.max()) if change.size else 0.0
    panels = [
        {
            "name": name,
            "data": [float(v) if np.isfinite(v) else None for v in a.ravel()],
            "low": lo,
            "high": hi,
        }
        for name, a, lo, hi in zip(
            ("Original terrain", "Healed terrain", "Applied difference"),
            arrays,
            (low, low, 0.0),
            (high, high, maximum),
        )
    ]
    payload = json.dumps(
        {"panels": panels, "width": int(window.width), "height": int(window.height)}
    )
    summary = html.escape(
        json.dumps(
            {
                "preview_window": {
                    "column": int(window.col_off),
                    "row": int(window.row_off),
                    "width": int(window.width),
                    "height": int(window.height),
                },
                "grid": qa["grid"],
                "vertical": qa["vertical"],
                "datum_check": qa["datum_check"],
                "counts": qa["counts"],
            },
            indent=2,
        )
    )
    document = """<!doctype html><html lang="en"><meta charset="utf-8">
<title>Mapping terrain review</title><style>
body{font:16px system-ui;margin:24px;background:#f4f5f7;color:#17232d}
main{display:flex;flex-wrap:wrap;gap:20px}canvas{width:384px;max-width:100%;image-rendering:pixelated;border:1px solid #777}
article{background:white;padding:12px}pre{white-space:pre-wrap}small{display:block}
</style><h1>Mapping terrain review</h1>
<p>Aligned pixels at native resolution. Original and healed share one elevation scale.
Difference shows applied raises. Grey pixels are NoData.</p>
<p>First modified block, up to 256 × 256 cells; if unchanged, the upper-left block.
Open the full TIFFs together in GIS to inspect other roads, bridges, channels and edges.
This preview is not engineering approval.</p><main></main><h2>QA record</h2><pre>SUMMARY</pre>
<script>const p=PAYLOAD;
for(const panel of p.panels){const article=document.createElement('article');
const title=document.createElement('h2');title.textContent=panel.name;article.append(title);
const c=document.createElement('canvas');c.width=p.width;c.height=p.height;article.append(c);
const ctx=c.getContext('2d'),im=ctx.createImageData(p.width,p.height);
panel.data.forEach((v,i)=>{const t=v===null?0:Math.max(0,Math.min(1,(v-panel.low)/(panel.high-panel.low||1)));
const rgb=v===null?[180,180,180]:[Math.round(240-190*t),Math.round(245-120*t),Math.round(250-80*t)];
im.data.set([...rgb,255],i*4)});ctx.putImageData(im,0,0);
const label=document.createElement('small');label.textContent=`Range: ${panel.low.toFixed(3)} to ${panel.high.toFixed(3)} model vertical units`;
article.append(label);document.querySelector('main').append(article);}
</script></html>""".replace("SUMMARY", summary).replace("PAYLOAD", payload)
    (output / "review.html").write_text(document, encoding="utf8")
