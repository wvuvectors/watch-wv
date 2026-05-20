CREATE DATABASE patchr_db;
USE patchr_db;

CREATE TABLE samples ( -- sample metadata
    sample_id VARCHAR(50) PRIMARY KEY,
    sample_status VARCHAR(50),
    location_id VARCHAR(100),
    sample_event VARCHAR(100),
    sample_qc VARCHAR(100),
    sample_collection_start_datetime DATETIME,
    sample_collection_end_datetime DATETIME,
    sample_recovered_datetime DATETIME,
    sample_collection_by VARCHAR(50),
    sample_flow FLOAT,
    sample_received_by VARCHAR(50),
    sample_received_date DATE,
    sample_ph_lab DECIMAL(4,3),
    sample_comment TEXT,
    FULLTEXT INDEX idx_comment (sample_comment),
    INDEX idx_location (location_id),
    INDEX idx_dates (sample_collection_start_datetime, sample_collection_end_datetime)
);

CREATE TABLE concentration ( -- concentration records
	concentration_id VARCHAR(50) PRIMARY KEY,
    sample_id VARCHAR(50),
    concentration_batch_id VARCHAR(50),
    concentration_location_in_batch VARCHAR(50),
    concentration_comment TEXT,
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id),
    INDEX idx_batch (concentration_batch_id),
    FULLTEXT INDEX idx_concentration_comment (concentration_comment)
);

CREATE TABLE cbatch ( -- concentration batch metadata
    concentration_batch_id VARCHAR(50) PRIMARY KEY,
    concentration_date DATE,
    concentration_input_ml FLOAT,
    concentration_machine VARCHAR(100),
    concentration_method VARCHAR(100),
    concentration_method_lot_id VARCHAR(100),
    concentration_output_ml FLOAT,
    concentration_run_by VARCHAR(100),
    concentration_batch_record_version VARCHAR(50),
    concentration_batch_comment TEXT,
    FULLTEXT INDEX idx_batch_comment (concentration_batch_comment),
    INDEX idx_concentration_date (concentration_date)
);

CREATE TABLE extractions ( -- individual extraction records
	extraction_id VARCHAR(50) PRIMARY KEY,
    concentration_id VARCHAR(50),
    extraction_batch_id VARCHAR(50),
    extraction_location_in_batch VARCHAR(50),
    extraction_location_in_storage VARCHAR(50),
    extraction_comment TEXT,
    sample_id VARCHAR(50),
    FOREIGN KEY (concentration_id) REFERENCES concentration(concentration_id),
    INDEX idx_batch (extraction_batch_id),
    FULLTEXT INDEX idx_extraction_comment (extraction_comment)
);

CREATE TABLE ebatch ( -- extraction batch metadata
	extraction_batch_id VARCHAR(50) PRIMARY KEY,
    extraction_date DATE,
    extraction_input_ul FLOAT,
    extraction_eluant VARCHAR(100),
    extraction_machine VARCHAR(100),
    extraction_method VARCHAR(100),
    extraction_method_lot_id VARCHAR(100),
    extraction_output_ul FLOAT,
    extraction_batch_record_version VARCHAR(50),
    extraction_run_by VARCHAR(100),
    extraction_batch_comment TEXT,
    FULLTEXT INDEX idx_batch_comment (extraction_batch_comment),
    INDEX idx_extraction_date (extraction_date)
);

CREATE TABLE assay ( -- assay run metadata
    assay_id VARCHAR(50) PRIMARY KEY,
    extraction_id VARCHAR(50),
    sample_id VARCHAR(50),
    assay_batch_id VARCHAR(50),
    assay_location_in_batch VARCHAR(50),
    assay_input_ul FLOAT,
    assay_class VARCHAR(100),
    assay_type VARCHAR(100),
    assay_target VARCHAR(200),
    assay_target_genetic_locus VARCHAR(100),
    assay_template VARCHAR(100),
    assay_target_macromolecule VARCHAR(100),
    assay_target_fluorophore VARCHAR(100),
    assay_accepted_droplets INT,
    assay_target_predicted_copies_per_ul_reaction DOUBLE,
    assay_target_copies_per_ul_reaction DOUBLE,
    assay_comment TEXT,
    FOREIGN KEY (extraction_id) REFERENCES extractions(extraction_id),
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id),
    INDEX idx_batch (assay_batch_id),
    INDEX idx_target (assay_target),
    FULLTEXT INDEX idx_assay_comment (assay_comment)
);

CREATE TABLE abatch ( -- assay batch metadata
    assay_batch_id VARCHAR(50) PRIMARY KEY,
    assay_date DATE,
    assay_reaction_ul FLOAT,
    assay_machine VARCHAR(100),
    assay_amplification_method VARCHAR(100),
    assay_amplification_method_lot_id VARCHAR(100),
    assay_quantification_method VARCHAR(100),
    assay_quantification_type VARCHAR(100),
    assay_batch_record_version VARCHAR(50),
    assay_method VARCHAR(255),
    assay_method_lot_id VARCHAR(100),
    assay_qx_manager_version VARCHAR(100),
    assay_run_by VARCHAR(100),
    assay_batch_comment TEXT,
    FULLTEXT INDEX idx_batch_comment (assay_batch_comment),
    INDEX idx_assay_date (assay_date)
);

CREATE TABLE results (
    assay_id VARCHAR(50) PRIMARY KEY,
    sample_id VARCHAR(50),
    collection_start_datetime DATETIME,
    collection_end_datetime DATETIME,
    event_type VARCHAR(100),
    sample_flow FLOAT,
    sample_qc VARCHAR(50),
    location_id VARCHAR(100),
    target VARCHAR(200),
    target_genetic_locus VARCHAR(100),
    lab_id VARCHAR(50),
    target_copies_per_l DOUBLE,
    target_copies_per_ld DOUBLE,
    target_copies_per_ldcap DOUBLE,
    target_copies_flownorm DOUBLE,
    target_copies_fn_per_cap DOUBLE,
    target_per_capita_basis DOUBLE,
    nc_copies_per_rxn DOUBLE,
    pc_copies_per_rxn DOUBLE,
    target_result_validated VARCHAR(100),
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id),
    FOREIGN KEY (assay_id) REFERENCES assay(assay_id),
    INDEX idx_sample (sample_id),
    INDEX idx_target (target),
    INDEX idx_location (location_id),
    FULLTEXT INDEX idx_validation (target_result_validated)
);

CREATE TABLE results_old (
    assay_id VARCHAR(50) PRIMARY KEY,
    sample_id VARCHAR(50),
    collection_start_datetime DATETIME,
    collection_end_datetime DATETIME,
    event_type VARCHAR(100),
    sample_flow FLOAT,
    sample_qc VARCHAR(50),
    location_id VARCHAR(100),
    target VARCHAR(200),
    target_genetic_locus VARCHAR(100),
    lab_id VARCHAR(50),
    target_copies_per_l DOUBLE,
    target_copies_per_ld DOUBLE,
    target_copies_per_ldcap DOUBLE,
    target_copies_flownorm DOUBLE,
    target_copies_fn_per_cap DOUBLE,
    target_per_capita_basis DOUBLE,
    nc_copies_per_rxn DOUBLE,
    pc_copies_per_rxn DOUBLE,
    target_result_validated VARCHAR(100),
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id),
    INDEX idx_sample (sample_id),
    INDEX idx_target (target),
    INDEX idx_location (location_id),
    FULLTEXT INDEX idx_validation (target_result_validated)
);

CREATE TABLE archive ( -- archive metadata for stored samples
    archive_id VARCHAR(50) PRIMARY KEY,
    sample_id VARCHAR(50),
    archive_batch_id VARCHAR(50),
    archive_location_in_batch VARCHAR(50),
    archive_location_in_storage VARCHAR(50),
    archive_comment TEXT,
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id),
    INDEX idx_batch (archive_batch_id),
    FULLTEXT INDEX idx_archive_comment (archive_comment)
);

CREATE TABLE rbatch ( -- archive batch metadata
    archive_batch_id VARCHAR(50) PRIMARY KEY,
    archive_date DATE,
    archive_input_ml FLOAT,
    archive_eluant VARCHAR(100),
    archive_machine VARCHAR(100),
    archive_method VARCHAR(100),
    archive_method_lot_id VARCHAR(100),
    archive_output_ml FLOAT,
    archive_run_by VARCHAR(100),
    archive_batch_record_version VARCHAR(50),
    archive_batch_comment TEXT,
    FULLTEXT INDEX idx_batch_comment (archive_batch_comment),
    INDEX idx_archive_date (archive_date)
);

CREATE TABLE files (
    file_id INT AUTO_INCREMENT PRIMARY KEY,
    filename VARCHAR(255),
    filepath VARCHAR(500),
    file_size INT,
    last_modified DATETIME,
    upload_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    INDEX idx_filename (filename)
);


-- TRUNCATE TABLE samples;
-- TRUNCATE TABLE extractions;
-- TRUNCATE TABLE ebatch;
-- TRUNCATE TABLE cbatch;
-- TRUNCATE TABLE concentration;
-- TRUNCATE TABLE abatch;
-- TRUNCATE TABLE rbatch;
-- TRUNCATE TABLE assay;
-- TRUNCATE TABLE results;
-- TRUNCATE TABLE result_old;
-- TRUNCATE TABLE archive;



    

