from pathlib import Path

#Having to fix relative paths all the time got tiring...

PROJECT_ROOT = Path(__file__).parent.parent.resolve()

INTERNAL_FILES_DIR = PROJECT_ROOT / "Internal_Files"
COMPUTATION_DIR = PROJECT_ROOT / "Computation_Deposit"
IMAGES_DIR = PROJECT_ROOT / "Images"
SQL_DIR = PROJECT_ROOT / "sql_files"