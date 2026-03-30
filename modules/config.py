import configparser

# Valid tasks
VALID_TASKS = ['readquality', 'filter', 'assembly', 'characterization', 'all']

# Config parameters: (key, type, default_value)
PARAM_SCHEMAS = {
    'filter': [
        ('cut_adapter', 'str', 'y'),
        ('filter_quality', 'int', 7),
        ('filter_length', 'int', 1000),
        ('porechop_threads', 'int', 8),
    ],
    'assembly': [
        ('flye_mode', 'str', 'nano-hq'),
        ('flye_meta', 'str', 'y'),
        ('flye_threads', 'int', 8),
        ('goldrush_genome_size', 'int', 300000),
        ('goldrush_min_read_length', 'int', 1000),
        ('goldrush_threads', 'int', 8),
        ('raven_threads', 'int', 8),
        ('shasta_read_length', 'int', 1000),
        ('shasta_threads', 'int', 8),
        ('wtdbg2_genome_size', 'int', 300000),
        ('wtdbg2_read_length', 'int', 1000),
        ('wtdbg2_kmer_fsize', 'int', 15),
        ('wtdbg2_kmer_psize', 'int', 0),
        ('wtdbg2_min_similarity', 'float', 0.05),
        ('wtdbg2_realign', 'str', 'y'),
        ('wtdbg2_threads', 'int', 8),
        ('polish_threads', 'int', 2),
        ('checkv_db', 'str', '~/dbs/checkv-db-v1.5'),
        ('checkv_threads', 'int', 8),
    ],
    'characterization': [
        ('pharokka_threads', 'int', 8),
        ('pharokka_db', 'str', '~/dbs/pharokka_db'),
        ('taxmyphage_threads', 'int', 8),
        ('taxmyphage_db', 'str', '~/dbs/taxmyphage_db'),
        ('phabox_db', 'str', '~/dbs/phabox_db_v2'),
        ('phabox_threads', 'int', 8),
        ('vhulk_location', 'str', '~/vHULK'),
        ('vhulk_threads', 'int', 8),
    ],
}

# Function to generate default config file from PARAM_SCHEMAS using default_value
def generate_default_config(config_file='qppl.conf'):
    lines = ['[General]', 'input_dir = input', 'output_dir = output', 'task = all', '']
    
    for task_name, schema in PARAM_SCHEMAS.items():
        lines.append(f'[{task_name.capitalize()}]')
        for key, _, default_value in schema:
            lines.append(f'{key} = {default_value}')
        lines.append('')
    
    with open(config_file, 'w') as file:
        file.write('\n'.join(lines))
    print(f"Default configuration file '{config_file}' generated.")

# Function to retrieve and convert config values based on type
def _get_config_value(config, section, key, value_type, default_value=None):
    try:
        getter = {'int': config.getint, 'float': config.getfloat}.get(value_type, config.get)
        value = getter(section, key, fallback=None)  # Get value without fallback first

        # If the value is None or an empty string, return default_value
        if value in (None, ''):
            return default_value

        # If the value is equal to the default_value, return default_value
        if value == default_value:
            return default_value

        return value  # Return the actual value if it's a valid value
    except (configparser.NoOptionError, configparser.NoSectionError):
        return default_value  # Return default_value if the key or section is missing
    except ValueError:
        return default_value  # Return default_value if conversion fails (empty string to int/float)

# Function to load task parameters based on schema, using default_value if key missing
def _load_task_params(config, schema, section):
    return {key: _get_config_value(config, section, key, value_type, default_value)
            for key, value_type, default_value in schema}

# Function to read and parse configuration file
def read_config(config_file):
    config = configparser.ConfigParser()
    config.read(config_file)

    try:
        input_dir = config.get('General', 'input_dir')
        output_dir = config.get('General', 'output_dir')
        task = config.get('General', 'task')

        if task not in VALID_TASKS:
            raise ValueError(f"Invalid task '{task}'. Accepted tasks are: {', '.join(VALID_TASKS)}")

        task_params = {task_name: _load_task_params(config, schema, task_name.capitalize())
                       for task_name, schema in PARAM_SCHEMAS.items()
                       if task in [task_name, 'all']}

        return input_dir, output_dir, task, task_params

    except (configparser.NoOptionError, configparser.NoSectionError) as e:
        raise ValueError(f"Missing configuration: {e}")