# HipHoshiningen — 16S rRNA Sequencing Report

## TODO
Split this into three separate package would be wise. 
1. Hip would be the base ( parsing, enrich data )
2. hoshi would be group base report generator.
3. ningen is invidiual based report generator.

## Quickstart

Generate a simple HTML report containing one table per TSV file:

```bash
hoshi html-report test_data/emu_outputs/*.tsv -o dist/emu-samples.html
```

Omit `-o` to print the HTML to stdout. Use `--title` to customise the page heading.
