# Hoshi — 16S rRNA Sequencing Report generator.

## TODO

## Quickstart

### Other...
### report-single
Generate a single HTML report for a single sample:

```bash
 hoshi report-single -n sample01 test_data/emu_output/test_01/sample01 -o dist/sample01_report.html
```

### report-multi
Generate an aggregate HTML report from multiple TSV files:

```bash
hoshi report-multi test_data/raw_out/kash_emu/*.tsv -o dist/emu-samples.html
```


