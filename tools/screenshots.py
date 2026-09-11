#!/usr/bin/env python3
"""Capture the screenshots used by the About / Documentation / Tutorial pages.

    Rscript -e 'shiny::runApp(".", port = 8932, host = "127.0.0.1")' &
    python3 tools/screenshots.py [--port 8932] [--out www/images]

Requires selenium and Firefox; geckodriver is fetched by Selenium Manager.
Re-run this after any change to the layout so the documentation never shows a
UI that no longer exists -- which is exactly what the inherited v1 PNGs did.

Each shot navigates by deep link (?page=&tool=), drives the real selectize
widgets so the sidebar and the figure agree, waits for plotly to finish, then
captures the full page.
"""
import argparse, os, sys, time

try:
    from selenium import webdriver
    from selenium.webdriver.firefox.options import Options
except ImportError:
    sys.exit("selenium is required:  python3 -m pip install --user selenium")

# --- the views to capture --------------------------------------------------
# js: statements run after the page settles, before the screenshot.
SEL = "$('#{id}')[0].selectize.setValue({val});"
CLICK = "$('#{id}').click();"

SHOTS = [
    dict(name="ui_input_upload", page="Analyse", tool="Input data",
         caption="The upload panel: prediction file, FASTA, detected format."),

    dict(name="ui_input_loaded", page="Run example", tool="Input data",
         caption="A loaded dataset: summary strip and the parsed peptide x allele table."),

    dict(name="tool_distribution", page="Run example", tool="Distribution",
         js=[SEL.format(id="demo_distribution-alleles",
                        val="['HLA-A01:01','HLA-A02:01','HLA-A03:01']"),
             CLICK.format(id="demo_distribution-go")],
         caption="Epitope Distribution: score histogram, epitope table, and the per-allele cutoff ladder."),

    dict(name="tool_intersection", page="Run example", tool="Intersection",
         js=[SEL.format(id="demo_intersection-alleles",
                        val="['HLA-A01:01','HLA-A02:01','HLA-A03:01','HLA-B07:02']"),
             CLICK.format(id="demo_intersection-go")],
         caption="Epitope Intersection: an UpSet plot of epitopes shared between MHC alleles."),

    dict(name="tool_intersection_venn", page="Run example", tool="Intersection",
         js=[SEL.format(id="demo_intersection-alleles",
                        val="['HLA-A01:01','HLA-A02:01','HLA-A03:01']"),
             CLICK.format(id="demo_intersection-go"),
             "$('input[name=\\'demo_intersection-plot_type\\'][value=\\'Venn diagram\\']')"
             ".prop('checked',true).trigger('change');"],
         caption="The same comparison as a Venn diagram, for two to four alleles."),

    dict(name="tool_density", page="Run example", tool="Density",
         caption="Epitope Density: count against protein length, plus the protein x allele grid."),

    dict(name="tool_viewer", page="Run example", tool="Viewer",
         js=[SEL.format(id="demo_viewer-protein", val="'sp|P0DTC5|VME1_SARS2'"),
             CLICK.format(id="demo_viewer-go")],
         caption="Epitope Viewer: epitopes placed along the membrane protein, coloured by alleles bound."),

    dict(name="tool_promiscuity", page="Run example", tool="Promiscuity",
         caption="Epitope Promiscuity: strong and weak binders for the broadest-binding epitopes."),

    dict(name="tool_conservation", page="Run example", tool="Conservation",
         js=[SEL.format(id="demo_conservation-proteins",
                        val="['sp|P0DTC2|SPIKE_SARS2','sp|P0DTC1|R1A_SARS2',"
                            "'sp|P0DTD1|R1AB_SARS2','sp|P0DTC9|NCAP_SARS2']"),
             CLICK.format(id="demo_conservation-go")],
         caption="Epitope Conservation: epitopes shared between proteins, strains or variants."),
]

READY = """
var n = document.querySelectorAll('.js-plotly-plot').length;
var busy = document.querySelectorAll('.recalculating, .shiny-busy').length;
return (n > 0 && busy === 0);
"""


def wait_ready(d, timeout=45):
    """Wait for plotly to have drawn and Shiny to be idle."""
    end = time.time() + timeout
    while time.time() < end:
        try:
            if d.execute_script(READY):
                time.sleep(1.2)          # let the final paint land
                return True
        except Exception:
            pass
        time.sleep(0.5)
    return False


def full_height(d, minimum=900):
    return max(minimum, d.execute_script(
        "return Math.max(document.body.scrollHeight, "
        "document.documentElement.scrollHeight);"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--port", type=int, default=8932)
    ap.add_argument("--out", default="www/images")
    ap.add_argument("--width", type=int, default=1600)
    ap.add_argument("--only", default=None, help="capture just this shot name")
    a = ap.parse_args()

    os.makedirs(a.out, exist_ok=True)
    base = f"http://127.0.0.1:{a.port}/"

    o = Options()
    o.add_argument("-headless")
    o.set_preference("layout.css.devPixelsPerPx", "1.0")
    d = webdriver.Firefox(options=o)
    d.set_window_size(a.width, 1200)

    try:
        for s in SHOTS:
            if a.only and s["name"] != a.only:
                continue
            url = f"{base}?page={s['page'].replace(' ', '+')}&tool={s['tool'].replace(' ', '+')}"
            d.get(url)
            if not wait_ready(d):
                print(f"  ! {s['name']}: timed out waiting for the page to settle")
            for js in s.get("js", []):
                d.execute_script(js)
                time.sleep(0.6)
            if s.get("js"):
                wait_ready(d)
            # Grow the viewport to the document so the whole page is in frame.
            d.set_window_size(a.width, full_height(d))
            time.sleep(1.0)
            path = os.path.join(a.out, s["name"] + ".png")
            d.save_screenshot(path)
            kb = os.path.getsize(path) // 1024
            print(f"  {s['name']:<24} {kb:>5} KB   {path}")
            d.set_window_size(a.width, 1200)
    finally:
        d.quit()


if __name__ == "__main__":
    main()
