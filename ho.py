from pathlib import Path
from jinja2 import Environment, FileSystemLoader

# Paths
BASE = Path(__file__).parent
TPL_DIR = BASE / "templates"

# Load Jinja2
env = Environment(loader=FileSystemLoader(str(TPL_DIR)), autoescape=False)
tpl = env.get_template("microbook.html")

# Sample context (edit as needed)
context = {
  "page_title":"16S rRNA Sequencing Report",
  "logo_text":"LOGO",
  "lab_header":"Long-read Lab Reference",
  "accession_id":"LRL-24-000123",
  "report_no":"16S-000987",
  "patient":{"name":"NA","hn":"H123456","dob":"NA","age":"NA","gender":"NA","ethnicity":"NA"},
  "dates":{"collected":"2025-02-21","received":"2025-02-22","tested":"2025-02-22","reported":"2025-02-24"},
  "requester":"NA", "provider":"Siriraj Hospital",
  "specimen":{"id":"S-01","type":"Blood (EDTA), suspected sepsis"},
  "test":{"required":"Pathogen detection (bacterial), culture-independent","performed":"16S rRNA amplicon sequencing (V1–V9), ONT"},
  "taxa":[
    {"name":"Streptococcus agalactiae","common_name":"Group B Streptococcus","reads":15842,"rel_abundance":"62.1%","confidence":"High","interpretation":"Significant—compatible with invasive infection","significant":True},
    {"name":"Gardnerella swidsinskii","reads":1214,"rel_abundance":"4.8%","confidence":"Moderate","interpretation":"Detected—clinical relevance uncertain outside GU sites","significant":False},
    {"name":"Gardnerella vaginalis","reads":6345,"rel_abundance":"24.9%","confidence":"High","interpretation":"Significant—consider genitourinary source","significant":True},
    {"name":"Cutibacterium acnes","reads":412,"rel_abundance":"1.6%","confidence":"Moderate","interpretation":"Common skin commensal; possible background/contaminant","significant":False}
  ],
  "other_low_level":{"threshold":"1%","reads":"~695","rel_abundance":"~2.6%"},
  "qc":[
    {"metric":"Total raw reads","observed":25860,"threshold":"\u2265 10,000","status":"PASS"},
    {"metric":"Mean read Q-score","observed":12.6,"threshold":"\u2265 10.0","status":"PASS"},
    {"metric":"Host (human) read fraction","observed":"3.2%","threshold":"\u2264 20%","status":"PASS"},
    {"metric":"Negative control contamination","observed":"Not detected","threshold":"None","status":"PASS"},
    {"metric":"Positive control (mock community)","observed":"Detected (expected profile)","threshold":"Detect all expected","status":"PASS"},
    {"metric":"Coverage of 16S variable regions","observed":"V1–V9","threshold":"Assay-defined","status":"PASS"}
  ],
  "reviewers":{"cls":"Miss Linda Clarke","second":"Mr. Mira Ashford","supervisor":"Miss Jordan Ember"},
  "recommendations":[
    "Correlate with blood culture/Gram stain and clinical context.",
    "Species-level calls within the Gardnerella genus may require targeted tests if clinically indicated.",
    "For suspected invasive disease, follow institutional antimicrobial guidelines."
  ],
  "technical":{"sequencer":"MinION (ONT)","kit":"16S Barcoding Kit","software":"nano16S v0.9; Kraken2/Bracken; SILVA 138.1"},
  "cover":{"title":"Microbook — 16S rRNA Sequencing Report","subtitle":"Bacterial Pathogen Detection from Clinical Specimens","footer":"Department of Medical Bioinformatics • Siriraj Hospital"},
  "appendix":{
    "processing":"DNA extraction with host depletion; controls processed alongside clinical specimens.",
    "library":"ONT 16S barcoding workflow; amplicon pooling by equal molarity.",
    "bioinfo":"Basecalling (HAC), trimming, human depletion, Kraken2/Bracken with SILVA 138.1; lab SOP thresholds."
  },
  "disclaimers":{"lab":"Long-read Lab Reference"}
}

# Render HTML
html = tpl.render(**context)
out_html = BASE / "microbook_rendered.html"
out_html.write_text(html, encoding="utf-8")
print(f"Wrote {out_html}")

# Optional: PDF generation with WeasyPrint (uncomment if you have it installed)
from weasyprint import HTML
HTML(filename=str(out_html)).write_pdf(str(BASE / "microbook.pdf"))
print("Wrote microbook.pdf")
