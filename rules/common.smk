
# exp_des = (pd.read_table(config["DE"]["Experimental_design"])
#              .set_index("Sample_ID", 
#                         drop=False))
# validate(exp_des, schema = "../schemas/Experimental_design.schema.yaml")

project_root = Path.cwd()

#### Shortened paths/names ####
# Results
FIGS = 'results/figs'
PICS = 'results/pics'
TBLS = 'results/tbls'

# Workflow paths
ENVS = '../envs'
#ENVS = project_root / 'envs' # Path() doesn't work https://github.com/snakemake/snakemake/issues/798#issue-762994335
REP = '../report/doc'
REP = project_root/'report'/'doc'
BMARKS = project_root/'bmarks'
SRC = project_root/'src'
LOGS = project_root/'logs'


# Data
RAW = 'data/raw'
REF = 'data/ref'
PROC = 'data/interm'

# Report categories
DE = "DE analysis results"
QC = "QC"
GSEA = "GSEA"
SUM = "Final results"

# Common variables for rule all targets
SHRINKAGE = ["not-shrinked","shrinked"]

# In case p-adj is given by the user convert it to float
padj_limit = float(config["DE"]["p-adj"])
