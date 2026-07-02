---
name: datafile-quicklook
description: Quickly inspect CSV/TSV/Excel and basic text data for schema, encoding, missing values, and anomalies. Use when the user mentions .csv/.tsv/.xls/.xlsx, data preview, or data cleaning.
---

# Datafile Quicklook

## What to extract
- File type, size, encoding (common: UTF-8/GBK), delimiter
- Columns: names, inferred types, missingness
- Basic stats: count, unique, min/max for numeric, top values for categorical
- Obvious issues: mixed types, whitespace, duplicated columns/rows, bad headers

## Safe workflow
1. Read a small sample first (head), then full file if needed.
2. Avoid loading huge Excel sheets fully; select sheet and columns.
3. If encoding is unclear, try UTF-8 then GBK.

## Output format
- **Schema**
- **Data quality findings**
- **Recommended cleaning steps**
