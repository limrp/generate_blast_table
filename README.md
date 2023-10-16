# generate_blast_table.py

## Description
`generate_blast_table.py` is a Python script designed to parse BLAST XML results and generate a summary table. This table is used to classify proteins into families based on the BLAST hits. The table is generated in Excel format, with an optional CSV output.

## Final Product
The script produces a DataFrame that classifies proteins by families. The DataFrame has the following columns:

- `Protein Family`: The name of the protein family.
- `Total BLAST Hits`: Total number of hits corresponding to each protein family.
- `Filtered BLAST Hits`: Number of hits that pass certain filtering criteria (e.g., E-value, alignment span).
- `Unique BLAST Hits`: Number of unique hits per protein family.

## Dependencies

The script requires the following Python packages:

- BioPython
- pandas
- alive-progress

You can install these dependencies using `pip`:

```bash
pip install biopython pandas alive-progress
```

## Usage

### Basic Usage

```bash
python generate_blast_table.py -i input.xml -o output.xlsx
```

### With All Options

```bash
python generate_blast_table.py -i input.xml -o output.xlsx -o2 output.csv -p custom_protein_families.txt -e 1e-5 -a 100 -l my_log.txt
```

### Arguments

- `-i, --input` (required): Input BLAST XML file.
- `-o, --output_excel` (required): Output Excel file.
- `-o2, --output_csv` (optional): Output CSV file.
- `-p, --protein_families` (optional): File containing a custom list of protein families. If not provided, the script will use a default list.
- `-e, --evalue` (optional): E-value threshold. Default is `1e-10`.
- `-a, --aln_span` (optional): Alignment span threshold. Default is `50`.
- `-l, --log` (optional): Log file name. Default is `log.txt`.


