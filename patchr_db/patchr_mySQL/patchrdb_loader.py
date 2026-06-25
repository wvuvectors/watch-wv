# patchrdb_loader.py
# Automatically loads tab-delimited watchdb files into MySQL,
# respecting foreign key order and skipping invalid child rows

import os
import mysql.connector
import pandas as pd

# ----------------------------- CONFIGURATION -----------------------------
DB_CONFIG = {
    'host': 'localhost',
    'user': 'root',
    'password': 'WVUvectors4209!!',  
    'database': 'patchr_db'
}

WATCHDB_DIR = os.path.dirname(os.path.abspath(__file__))
LOG_FILE = os.path.join(WATCHDB_DIR, 'patchrdb_loader.log')
INVALID_DIR = os.path.join(WATCHDB_DIR, 'invalid_rows')

# -----------------------------
# Table mapping and order
# -----------------------------
table_map = {
    'sample': 'samples',
    'concentration': 'concentration',
    'extraction': 'extractions',
    'assay': 'assay',
    'result.txt': 'results',
    'result.OLD': 'results_old',
    'cbatch': 'cbatch',
    'ebatch': 'ebatch',
    'abatch': 'abatch',
    'rbatch': 'rbatch',
    'archive': 'archive',
    'location': 'location',
}

table_order = [
    'samples',
    'cbatch',
    'concentration',
    'ebatch',
    'extractions',
    'abatch',
    'assay',
    'results',
    'results_old',
    'archive',
    'rbatch',
    'location'
]

# -----------------------------
# Database connection
# -----------------------------
def connect_db():
    return mysql.connector.connect(**DB_CONFIG)

def log_message(msg):
    print(msg)
    with open(LOG_FILE, 'a', encoding='utf-8') as f:
        f.write(msg + '\n')

# -----------------------------
# Data cleaning
# -----------------------------
import numpy as np

def clean_dataframe(df):
    # --------------------------------------------------
    # TEMPORARY FIX:
    # watchdb.abatch.txt renamed assay_qx_manager_version
    # to assay_analysis_software_version
    # --------------------------------------------------
    if (
        "assay_analysis_software_version" in df.columns
        and "assay_qx_manager_version" not in df.columns
    ):
        log_message(
            "[INFO] Renaming assay_analysis_software_version "
            "-> assay_qx_manager_version for compatibility"
        )

        df = df.rename(
            columns={
                "assay_analysis_software_version":
                "assay_qx_manager_version"
            }
        )

    
    # Clean ID, date/datetime, and numeric columns before insertion.
    for col in df.columns:
        col_lower = col.lower()
        
        # ID columns: strip & uppercase
        if 'id' in col_lower:
            df[col] = df[col].astype(str).str.strip().str.upper()
            df[col] = df[col].replace({'NAN': None, '': None})
        
        # Date columns: parse to YYYY-MM-DD
        elif 'datetime' in col_lower:
            df[col] = pd.to_datetime(df[col], errors='coerce', format='mixed').dt.strftime('%Y-%m-%d %H:%M:%S')
            failed_rows = df[col].isna().sum()
            if failed_rows > 0: 
                log_message(f"[WARN] {failed_rows} values in "
                            f"{col} could not be parsed as datetime")
        #elif 'datetime' in col.lower():
            # Keep time information
        #    df[col] = pd.to_datetime(df[col], errors='coerce').dt.strftime('%Y-%m-%d %H:%M:%S')
        elif 'date' in col.lower():
            # Only keep the date portion
            df[col] = pd.to_datetime(df[col], errors='coerce').dt.strftime('%Y-%m-%d')
        
        # Numeric columns: convert to float/decimal
        elif 'ph' in col_lower or 'flow' in col_lower or 'copies' in col_lower or 'ml' in col_lower or 'ul' in col_lower:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    
    # Replace all remaining NaN with None for MySQL
    df = df.where(pd.notna(df), None)
    
    return df

def export_invalid_rows(df_invalid, table_name, reason):
    """Save invalid rows to CSV for review."""
    if df_invalid.empty:
        return
    os.makedirs(INVALID_DIR, exist_ok=True)
    filename = f"invalid_{table_name}_{reason}.csv"
    out_path = os.path.join(INVALID_DIR, filename)
    try:
        df_invalid.to_csv(out_path, index=False)
        log_message(f"[EXPORT] Saved {len(df_invalid)} invalid {table_name} rows to {filename}")
    except PermissionError as e:
        log_message(f"[ERROR] Could not export invalid rows: {e}")

# -----------------------------
# Foreign key validation
# -----------------------------
def validate_foreign_keys(df, table_name, cursor, export_invalid=True):
    """
    Remove rows whose parent IDs don't exist.
    - df: DataFrame to validate
    - table_name: current table
    - cursor: MySQL cursor
    - export_invalid: if True, save invalid rows to CSV
    """
    fk_checks = {
        'concentration': ('sample_id', 'samples', 'sample_id'),
        'extractions': ('concentration_id', 'concentration', 'concentration_id'),
        'assay': ('extraction_id', 'extractions', 'extraction_id'),
        'results': ('assay_id', 'assay', 'assay_id'),
        'results_old': ('sample_id', 'samples', 'sample_id'),
        'archive': ('sample_id', 'samples', 'sample_id'),
        'location': ('location_id', 'samples', 'location_id')
    }

    if table_name not in fk_checks:
        return df

    child_col, parent_table, parent_col = fk_checks[table_name]

    if child_col not in df.columns:
        return df

    # Fetch valid parent IDs
    cursor.execute(f"SELECT {parent_col} FROM {parent_table}")
    valid_ids = {row[0] for row in cursor.fetchall()}

    # Treat None/NaN as invalid
    df_valid = df[df[child_col].notna() & df[child_col].isin(valid_ids)]
    df_invalid = df[~df.index.isin(df_valid.index)]

    if not df_invalid.empty:
        log_message(f"[WARN] Skipping {len(df_invalid)} invalid {table_name} rows (missing {child_col}).")
        if export_invalid:
            export_invalid_rows(df_invalid, table_name, f'missing_{child_col}')

    return df_valid

# -----------------------------
# Insert function
# -----------------------------
def insert_dataframe(df, table_name, cursor):
    """
    Insert dataframe into MySQL with ON DUPLICATE KEY UPDATE.
    Handles None values correctly (instead of string 'nan').
    """
    if df.empty:
        return

    placeholders = ', '.join(['%s'] * len(df.columns))
    columns = ', '.join(df.columns)
    update_clause = ', '.join([f"{col}=VALUES({col})" for col in df.columns])
    sql = f"INSERT INTO {table_name} ({columns}) VALUES ({placeholders}) ON DUPLICATE KEY UPDATE {update_clause}"

    for _, row in df.iterrows():
        # Convert NaN or None to proper Python None
        values = tuple(None if pd.isna(v) else v for v in row.values)
        cursor.execute(sql, values)

# -----------------------------
# Main loader
# -----------------------------
def main():
    conn = connect_db()
    cursor = conn.cursor()
    log_message("\n========== PATCHR DB LOADER RUN ==========")
    print("WATCHDB_DIR =", WATCHDB_DIR)

    for table_name in table_order:
        file_to_load = None
        for file in os.listdir(WATCHDB_DIR):
            if not file.startswith('watchdb.') or not file.endswith('.txt'):
                continue
            for key, tbl in table_map.items():
                if key in file and tbl == table_name:
                    file_to_load = file
                    break
            if file_to_load:
                break

        if not file_to_load:
            log_message(f"[SKIP] No file found for table {table_name}")
            continue

        file_path = os.path.join(WATCHDB_DIR, file_to_load)
        log_message(f"[LOAD] {file_to_load} → {table_name}")

        try:
            df = pd.read_csv(file_path, sep='\t', dtype=str)
            df = clean_dataframe(df)
            df = validate_foreign_keys(df, table_name, cursor)

            if df.empty:
                log_message(f"[WARN] No valid rows left for {table_name}, skipping.")
                continue

            insert_dataframe(df, table_name, cursor)
            conn.commit()
            log_message(f"[OK] Inserted/Updated {len(df)} rows into {table_name}")

        except Exception as e:
            conn.rollback()
            log_message(f"[ERROR] Failed to insert {file_to_load} → {e}")

    cursor.close()
    conn.close()
    log_message("========== LOADER COMPLETE ==========\n")

if __name__ == '__main__':
    main()
