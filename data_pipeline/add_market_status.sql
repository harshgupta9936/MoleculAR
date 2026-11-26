-- Migration script to add 'market_status' to the drugs table
-- This replaces a flat purge/exclude strategy with descriptive tagging.

ALTER TABLE drugs
ADD COLUMN market_status VARCHAR(50) DEFAULT 'Active US';

-- Example usage (to be populated by categorization script output):
-- UPDATE drugs SET market_status = 'Veterinary' WHERE chembl_id = '...';
-- UPDATE drugs SET market_status = 'Withdrawn' WHERE chembl_id = '...';
-- UPDATE drugs SET market_status = 'Naming Mismatch' WHERE chembl_id = '...';
-- UPDATE drugs SET market_status = 'International/Non-US' WHERE chembl_id = '...';
