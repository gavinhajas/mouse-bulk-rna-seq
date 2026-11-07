#!/usr/bin/env python3
"""
NCBI & SRA Metadata Fetcher
Searches NCBI Entrez for SRA run accessions (SRR/ERR/DRR) 
or general NCBI accession numbers (RefSeq/GenBank) and prints key metadata.
"""

import sys
import xml.etree.ElementTree as ET
from Bio import Entrez

# ALWAYS set your email address for NCBI Entrez API guidelines
Entrez.email = "gavinsenpai1983@gmail.com"  # Replace with your email address
# Entrez.api_key = "YOUR_NCBI_API_KEY"  # Optional: set if you have an NCBI API key


def fetch_sra_info(sra_id: str):
    """Fetch metadata for an SRA accession number."""
    print(f"[*] Querying SRA database for: {sra_id} ...\n")
    try:
        # Search SRA database for ID
        handle = Entrez.esearch(db="sra", term=sra_id)
        search_results = Entrez.read(handle)
        handle.close()

        id_list = search_results.get("IdList", [])
        if not id_list:
            print(f"[!] No SRA record found for '{sra_id}'.")
            return

        # Fetch XML summary for SRA record
        handle = Entrez.efetch(db="sra", id=id_list[0], rettype="full", retmode="xml")
        xml_data = handle.read()
        handle.close()

        # Parse XML
        root = ET.fromstring(xml_data)
        
        # Extract details
        run_element = root.find(".//RUN")
        sample_element = root.find(".//SAMPLE")
        organism_element = root.find(".//ORGANISM")
        library_element = root.find(".//LIBRARY_DESCRIPTOR")
        platform_element = root.find(".//PLATFORM")

        title = root.findtext(".//STUDY_TITLE") or root.findtext(".//TITLE") or "N/A"
        organism = organism_element.findtext("SCIENTIFIC_NAME") if organism_element is not None else "N/A"
        strategy = library_element.findtext("LIBRARY_STRATEGY") if library_element is not None else "N/A"
        source = library_element.findtext("LIBRARY_SOURCE") if library_element is not None else "N/A"
        selection = library_element.findtext("LIBRARY_SELECTION") if library_element is not None else "N/A"
        
        # Extract Platform details
        platform_name = "N/A"
        model_name = "N/A"
        if platform_element is not None and len(platform_element) > 0:
            platform_name = platform_element[0].tag
            model_name = platform_element[0].findtext("INSTRUMENT_MODEL") or "N/A"

        # Read counts and spot sizes
        spots = run_element.attrib.get("total_spots", "N/A") if run_element is not None else "N/A"
        bases = run_element.attrib.get("total_bases", "N/A") if run_element is not None else "N/A"

        print("=" * 60)
        print(f" SRA METADATA: {sra_id}")
        print("=" * 60)
        print(f"Title:            {title.strip()}")
        print(f"Organism:         {organism}")
        print(f"Library Strategy: {strategy}")
        print(f"Library Source:   {source}")
        print(f"Library Selection:{selection}")
        print(f"Platform:         {platform_name} ({model_name})")
        print(f"Total Spots/Reads:{spots}")
        print(f"Total Bases:      {bases}")
        print("=" * 60)

    except Exception as e:
        print(f"[!] Error fetching SRA metadata: {e}")


def fetch_general_accession_info(accession: str):
    """Fetch metadata for a standard NCBI accession (GenBank/RefSeq/Nucleotide/Protein)."""
    print(f"[*] Querying NCBI Nucleotide/Protein database for: {accession} ...\n")
    try:
        # Search Nucleotide database first
        handle = Entrez.esearch(db="nuccore", term=accession)
        search_results = Entrez.read(handle)
        handle.close()

        db_used = "nuccore"
        id_list = search_results.get("IdList", [])

        # Fallback to Protein if not found in Nucleotide
        if not id_list:
            handle = Entrez.esearch(db="protein", term=accession)
            search_results = Entrez.read(handle)
            handle.close()
            id_list = search_results.get("IdList", [])
            db_used = "protein"

        if not id_list:
            print(f"[!] Accession '{accession}' not found in Nucleotide or Protein databases.")
            return

        # Fetch Summary
        handle = Entrez.esummary(db=db_used, id=id_list[0])
        summary = Entrez.read(handle)[0]
        handle.close()

        print("=" * 60)
        print(f" NCBI ACCESSION METADATA: {accession} ({db_used.upper()})")
        print("=" * 60)
        print(f"Accession:   {summary.get('AccessionVersion', accession)}")
        print(f"Title:       {summary.get('Title', 'N/A')}")
        print(f"Organism:    {summary.get('TaxId', 'N/A')} (TaxID)")
        print(f"Length:      {summary.get('Length', 'N/A')} bp/aa")
        print(f"Create Date: {summary.get('CreateDate', 'N/A')}")
        print(f"Update Date: {summary.get('UpdateDate', 'N/A')}")
        print("=" * 60)

    except Exception as e:
        print(f"[!] Error fetching accession metadata: {e}")


def lookup_id(query_id: str):
    """Detect whether query is an SRA run or general accession and route request."""
    query_clean = query_id.strip().upper()
    
    # SRA runs typically start with SRR, ERR, DRR, SRX, ERX, DRX, SRP, ERP, DRP
    if any(query_clean.startswith(prefix) for prefix in ["SRR", "ERR", "DRR", "SRX", "ERX", "DRX", "SRP", "ERP", "DRP"]):
        fetch_sra_info(query_clean)
    else:
        fetch_general_accession_info(query_clean)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python ncbi_lookup.py <ACCESSION_OR_SRA_ID>")
        print("Examples:")
        print("  python ncbi_lookup.py SRR1553610")
        print("  python ncbi_lookup.py NC_000913")
        sys.exit(1)

    user_query = sys.argv[1]
    lookup_id(user_query)
