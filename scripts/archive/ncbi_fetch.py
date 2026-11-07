#!/usr/bin/env python3
import json
import sys
import urllib.parse
import urllib.request


def print_progress_bar(iteration: int, total: int, length: int = 30):
    """Draws a visual text progress bar in the terminal."""
    percent = (iteration / total) * 100
    filled_length = int(length * iteration // total)
    bar = "█" * filled_length + "░" * (length - filled_length)
    sys.stdout.write(f"\r[{bar}] {percent:5.1f}% ({iteration}/{total})")
    sys.stdout.flush()


def search_ncbi(db: str, query: str):
    """Searches NCBI, uses WebEnv history, and returns count, WebEnv, and QueryKey."""
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": db,
        "term": query,
        "usehistory": "y",
        "retmode": "json",
    }
    url = f"{base_url}?{urllib.parse.urlencode(params)}"

    req = urllib.request.Request(url, headers={"User-Agent": "NCBIFetchScript/1.0"})
    with urllib.request.urlopen(req) as resp:
        data = json.loads(resp.read().decode("utf-8"))

    result = data.get("esearchresult", {})
    count = int(result.get("count", 0))
    webenv = result.get("webenv")
    query_key = result.get("querykey")

    return count, webenv, query_key


def fetch_and_save(
    db: str,
    rettype: str,
    webenv: str,
    query_key: str,
    count: int,
    output_file: str,
    batch_size: int = 100,
):
    """Fetches records in batches and displays a progress bar."""
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    with open(output_file, "w", encoding="utf-8") as out:
        # Initial empty bar
        print_progress_bar(0, count)

        for start in range(0, count, batch_size):
            # Calculate size for this specific slice
            current_batch = min(batch_size, count - start)

            params = {
                "db": db,
                "query_key": query_key,
                "WebEnv": webenv,
                "rettype": rettype,
                "retmode": "text",
                "retstart": start,
                "retmax": current_batch,
            }
            url = f"{base_url}?{urllib.parse.urlencode(params)}"
            req = urllib.request.Request(
                url, headers={"User-Agent": "NCBIFetchScript/1.0"}
            )

            with urllib.request.urlopen(req) as resp:
                chunk = resp.read().decode("utf-8")
                # Remove header lines on subsequent batches for CSV files
                if rettype == "runinfo" and start > 0:
                    lines = chunk.splitlines(keepends=True)
                    if lines:
                        chunk = "".join(lines[1:])
                out.write(chunk)

            # Update progress bar
            print_progress_bar(start + current_batch, count)

    print(f"\nDone! File saved to: {output_file}")


def main():
    query = input("Enter search query (e.g. PRJNA257197): ").strip()
    if not query:
        print("Error: Search query cannot be empty.")
        return

    print("\nChoose database:")
    print("  1: Nucleotide / DNA / RNA (FASTA)")
    print("  2: Protein / Amino Acid   (FASTA)")
    print("  3: SRA RunInfo            (CSV)")
    db_choice = (
        input("Select [1, 2, or 3] (default: 1): ").strip() or "1"
    )

    if db_choice == "2":
        db, rettype, output_path = "protein", "fasta", "proteins.fa"
    elif db_choice == "3":
        db, rettype, output_path = "sra", "runinfo", "runinfo.csv"
    else:
        db, rettype, output_path = "nucleotide", "fasta", "genomes.fa"

    custom_output = input(f"Enter output filename (press Enter for default '{output_path}'): ").strip()
    if custom_output:
        output_path = custom_output

    print(f"\nSearching {db} for '{query}'...")
    total_found, webenv, query_key = search_ncbi(db, query)

    if total_found == 0:
        print("No records found.")
        return

    print(f"Found {total_found} total record(s).")

    # Ask for record limit (defaults to 300)
    limit_input = input("How many to download? (Press Enter for 300, or type 'all'): ").strip()
    if limit_input.lower() == "all":
        download_count = total_found
    elif limit_input.isdigit():
        download_count = min(total_found, int(limit_input))
    else:
        download_count = min(total_found, 300)

    print(f"Downloading {download_count} record(s)...\n")
    # Batch size set to 100 so the progress bar updates incrementally
    fetch_and_save(db, rettype, webenv, query_key, download_count, output_path, batch_size=100)


if __name__ == "__main__":
    main()
