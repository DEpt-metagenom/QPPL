import configparser

VALID_TASKS = ['readquality', 'filter', 'assembly', 'characterization','all']

def generate_default_config(config_file='qppl.conf'):
    default_config = """[General]
input_dir = input
output_dir = output
prefix = qppl_run
task = all

[ReadQuality]

[Filter]
porechop_cut_adapter = y
nanotools_filter_quality = 7
nanotools_filter_length = 1000

[Assembly]
flye_mode = nano-hq
flye_meta = y
flye_threads = 8
goldrush_genome_size = 300000
goldrush_min_read_length = 1000
goldrush_threads = 8
raven_threads = 8
shasta_threads = 8
wtdbg2_genome_size = 300000
wtdbg2_read_length = 1000
wtdbg2_kmer_fsize = 15
wtdbg2_kmer_psize = 0
wtdbg2_min_similarity = 0.05
wtdbg2_realign = y
wtdbg2_threads = 8
checkv_db = ~/dbs/checkv-db-v1.5
checkv_threads = 8

[Characterization]
pharokka_threads = 8
pharokka_db = ~/dbs/pharokka_db
taxmyphage_threads = 8
phabox_db = ~/dbs/phabox_db_v2
phabox_threads = 8
vhulk_location = ~/vHULK
vhulk_threads = 8
"""
    with open(config_file, 'w') as file:
        file.write(default_config)
    
    print(f"Default configuration file '{config_file}' generated.")

def read_config(config_file):
    config = configparser.ConfigParser()
    config.read(config_file)

    try:
        # General Section
        input_dir = config.get('General', 'input_dir')
        output_dir = config.get('General', 'output_dir')
        prefix = config.get('General', 'prefix')
        task = config.get('General', 'task')

        if task not in VALID_TASKS:
            raise ValueError(f"Invalid task '{task}'. Accepted tasks are: {', '.join(VALID_TASKS)}")

        task_params = {}

        if task in ['filter', 'all']:
            # Filter Section parameters
            filter_params = {
                'cut_adapter': config.get('Filter', 'porechop_cut_adapter', fallback='y'),
                'filter_quality': config.getint('Filter', 'nanotools_filter_quality'),
                'filter_length': config.getint('Filter', 'nanotools_filter_length')
            }
            task_params['filter'] = filter_params

        if task in ['assembly', 'all']:
            # Assembly Section parameters
            assembly_params = {
                'flye_mode': config.get('Assembly', 'flye_mode'),
                'flye_meta': config.get('Assembly', 'flye_meta', fallback='y'),
                'flye_threads': config.getint('Assembly', 'flye_threads', fallback=1),
                'goldrush_genome_size': config.getint('Assembly', 'goldrush_genome_size'),
                'goldrush_min_read_length': config.getint('Assembly', 'goldrush_min_read_length'),
                'goldrush_threads': config.getint('Assembly', 'goldrush_threads', fallback=1),
                'raven_threads': config.getint('Assembly', 'raven_threads', fallback=1),
                'shasta_threads': config.getint('Assembly', 'shasta_threads', fallback=1),
                'wtdbg2_genome_size': config.getint('Assembly', 'wtdbg2_genome_size'),
                'wtdbg2_read_length': config.getint('Assembly', 'wtdbg2_read_length'),
                'wtdbg2_kmer_fsize': config.getint('Assembly', 'wtdbg2_kmer_fsize'),
                'wtdbg2_kmer_psize': config.getint('Assembly', 'wtdbg2_kmer_psize'),
                'wtdbg2_min_similarity': config.getfloat('Assembly', 'wtdbg2_min_similarity'),
                'wtdbg2_realign': config.get('Assembly', 'wtdbg2_realign', fallback='y'),
                'wtdbg2_threads': config.getint('Assembly', 'wtdbg2_threads', fallback=1),
                'checkv_db': config.get('Assembly', 'checkv_db'),
                'checkv_threads': config.getint('Assembly', 'checkv_threads', fallback=1)
            }
            task_params['assembly'] = assembly_params

        if task in ['characterization', 'all']:
            # Characterization Section parameters
            characterization_params = {
                'pharokka_threads': config.getint('Characterization', 'pharokka_threads', fallback=1),
                'pharokka_db': config.get('Characterization', 'pharokka_db'),
                'taxmyphage_threads': config.getint('Characterization', 'taxmyphage_threads', fallback=1),
                'phabox_db': config.get('Characterization', 'phabox_db'),
                'phabox_threads': config.getint('Characterization', 'phabox_threads', fallback=1),
                'vhulk_location': config.get('Characterization', 'vhulk_location'),
                'vhulk_threads': config.getint('Characterization', 'vhulk_threads', fallback=1)
            }
            task_params['characterization'] = characterization_params

        return input_dir, output_dir, prefix, task, task_params

    except (configparser.NoOptionError, configparser.NoSectionError) as e:
        raise ValueError(f"Missing configuration: {e}")
