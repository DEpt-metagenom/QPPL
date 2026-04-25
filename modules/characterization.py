import os
import subprocess
import re
import logging
import time

# Call logger
logger = logging.getLogger("QPPL")

# Output subdirectories for CHARACTERIZATION
OUTPUT_CHAR_DIRS = {
    1: "out_pharokka",
    2: "out_taxmyphage",
    3: "out_phabox",
    4: "out_vhulk",
    5: "out_abricate",
    6: "intermediate"
}

# Output subdirectories RESULTS
OUTPUT_RES_DIRS = {
    1: "genomes"
}

# All assemblers used
ASSEMBLERS = {
    1: "flye",
    2: "goldrush",
    3: "raven",
    4: "shasta",
    5: "wtdbg2"
}

# Function to fully log messages only to file
def log_to_file_only(message, level=logging.INFO):
    for h in logger.handlers:
        if isinstance(h, logging.FileHandler):
            record = logger.makeRecord(logger.name, level, None, None, message, None, None)
            h.emit(record)

# Function to set up logging
def setup_logging(log_path):
    file_handler = logging.FileHandler(log_path)
    file_handler.setFormatter(logging.Formatter("[%(levelname)s] %(asctime)s - %(message)s", "%Y-%m-%d %H:%M:%S"))

    # Remove previous file handlers to avoid duplicate logs
    for h in logger.handlers[:]:
        if isinstance(h, logging.FileHandler):
            logger.removeHandler(h)
    logger.addHandler(file_handler)

# Function to run a command in a conda environment
def run_command(env_name, command, file, env_number):
    log_to_file_only(f"Running command in {env_name} on {file}: $ {command}", logging.INFO)
    result = subprocess.run(["conda", "run", "-n", env_name, "bash", "-c", command], text=True, capture_output=True)

    if result.returncode != 0:
        logger.error(f"Tool failed for file: {file} (step {env_number})\nCommand: {command}\nError: {result.stderr.strip()}")
        return
    
    out = result.stdout.strip()
    if not out:
        return

    is_version_line = "Tool version" in command
    log_fn = logger.debug if is_version_line else logger.info
    for line in out.splitlines():
        log_fn(line)

# Function to generate commands for assembly tasks per environment
def generate_commands(characterization_dir, output_dir, file, task_params, basename, step_number):

    assembly_dir = os.path.join(output_dir, "assembly_task")

    commands_by_step = {
        1: [  # Pharokka
            f"echo \"Tool version (pharokka):\" $(pharokka.py --version)",
            rf"""if [ $(grep -c ">" {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes.fasta) -gt 1 ]; then pharokka.py -i {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes.fasta -o {characterization_dir}/out_pharokka/{basename} -t {task_params['characterization']['pharokka_threads']} -m --dnaapler -p {basename} -d {task_params['characterization']['pharokka_db']} -g prodigal-gv; else pharokka.py -i {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes.fasta -o {characterization_dir}/out_pharokka/{basename} -t {task_params['characterization']['pharokka_threads']} --dnaapler -p {basename} -d {task_params['characterization']['pharokka_db']} -g prodigal-gv; fi""",
            rf"""awk '/^LOCUS/ {{ if (out) close(out); out = "{characterization_dir}/out_pharokka/{basename}/" substr($2, 1) ".gbk"; print "Creating file:", out; }} {{ print > out }}' {characterization_dir}/out_pharokka/{basename}/{basename}.gbk"""
        ],
        2: [ # taxmyPhage
            f"echo \"Tool version (taxmyphage):\" $(conda list | awk '$1==\"taxmyphage\" {{print $2}}')",
            rf"taxmyphage run -i {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes.fasta -t {task_params['characterization']['taxmyphage_threads']} -db {task_params['characterization']['taxmyphage_db']} -o {characterization_dir}/{OUTPUT_CHAR_DIRS[2]}/{basename}"
        ],
        3: [ # PhaBOX2
            f"echo \"Tool version (phabox2):\" $(conda list | awk '$1==\"phabox\" {{print $2}}')",
            rf"phabox2 --task end_to_end --contigs {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes.fasta --outpth {characterization_dir}/{OUTPUT_CHAR_DIRS[3]}/{basename} --dbdir {task_params['characterization']['phabox_db']} --threads {task_params['characterization']['phabox_threads']}"
        ],
        4: [ # vHULK
            f"echo \"Tool version (vHULK):\" $(cd ~/vHULK && python ~/vHULK/vHULK.py --version && cd - >/dev/null)",
            rf"mkdir {task_params['characterization']['vhulk_location']}/{basename}",
            rf"mkdir {task_params['characterization']['vhulk_location']}/{basename}_in",
            rf"seqkit split -i {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes.fasta -O {task_params['characterization']['vhulk_location']}/{basename}_in --extension '.fasta'",
            rf"cd {task_params['characterization']['vhulk_location']} && python vHULK.py -i {basename}_in -o {basename} -t {task_params['characterization']['vhulk_threads']} && cd -",
            rf"mv {task_params['characterization']['vhulk_location']}/{basename} {characterization_dir}/out_vhulk/{basename}",
            rf"rm -r {task_params['characterization']['vhulk_location']}/{basename}_in",
            rf"echo 'contig,pred_genus,score_genus,entropy_genus,pred_species,score_species,entropy_species' > {characterization_dir}/out_vhulk/{basename}/{basename}_vHULK.csv",
            rf"""
            contigs=$(ls {characterization_dir}/out_vhulk/{basename}/predictions/* | sed -n 's/.*part_\([^\.]*\)\.csv/\1/p') &&
            for contig in $contigs; do
                tail -n +2 {characterization_dir}/out_vhulk/{basename}/predictions/*${{contig}}*.csv |
                awk -v sample=$contig '{{print sample "," $0}}' >> {characterization_dir}/out_vhulk/{basename}/{basename}_vHULK.csv
            done
            """
        ],
        5: [ # ABRICATE
            f"echo \"Tool version (abricate):\" $(abricate --version)",
            rf"mkdir {characterization_dir}/out_abricate/{basename}",
            rf"abricate -db vfdb {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes.fasta --quiet > {characterization_dir}/out_abricate/{basename}/vfdb_{basename}.abr",
            rf"abricate -db card {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes.fasta --quiet > {characterization_dir}/out_abricate/{basename}/card_{basename}.abr",
            rf"abricate -db ncbi {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes.fasta --quiet > {characterization_dir}/out_abricate/{basename}/ncbi_{basename}.abr",
            rf"abricate -db megares {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes.fasta --quiet > {characterization_dir}/out_abricate/{basename}/megares_{basename}.abr",
            rf"abricate --summary {characterization_dir}/out_abricate/{basename}/vfdb_{basename}.abr > {characterization_dir}/out_abricate/{basename}/vfdb_{basename}.tab",
            rf"abricate --summary {characterization_dir}/out_abricate/{basename}/card_{basename}.abr > {characterization_dir}/out_abricate/{basename}/card_{basename}.tab",
            rf"abricate --summary {characterization_dir}/out_abricate/{basename}/ncbi_{basename}.abr > {characterization_dir}/out_abricate/{basename}/ncbi_{basename}.tab",
            rf"abricate --summary {characterization_dir}/out_abricate/{basename}/megares_{basename}.abr > {characterization_dir}/out_abricate/{basename}/megares_{basename}.tab"
        ],
        6: [ # INTERMEDIATE FILES
            rf"""awk 'BEGIN {{FS=OFS=","}} NR==1 {{$1="contig"; $2="vclust_cluster_id"; $NF="repres_score"; for (i=3; i<=NF-1; i++) $i="checkv_"$i}} 1' {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes.csv | sed 's/,/\t/g' > {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes_final.tsv""",
            #
            rf"""sort -t $'\t' -k1,1 {characterization_dir}/out_pharokka/{basename}/{basename}_length_gc_cds_density.tsv | awk 'BEGIN {{FS=OFS="\t"}} NR==1 {{$1="contig"; for (i=2; i<=NF; i++) $i="pharokka_"$i}} 1' > {characterization_dir}/out_pharokka/{basename}/{basename}_length_gc_cds_density_sorted.tsv""",
            rf"""sort -t $'\t' -k1,1 {characterization_dir}/out_pharokka/{basename}/{basename}_top_hits_mash_inphared.tsv | awk 'BEGIN {{FS=OFS="\t"}} NR==1 {{$1="contig"; for (i=2; i<=NF; i++) $i="pharokka_closest_"$i}} 1' > {characterization_dir}/out_pharokka/{basename}/{basename}_top_hits_mash_inphared_sorted.tsv""",
            rf"""{{
            printf "contig\tpharokka_CDS\tpharokka_tRNAs\tpharokka_tmRNAs\tpharokka_CRISPRs\n";
            awk 'BEGIN {{FS=OFS="\t"}}
                NR==1 {{next}}
                {{
                    d=$1; c=$2+0; ctg=$3;
                    if (d=="CDS") cds[ctg]=c;
                    else if (d=="tRNAs") trnas[ctg]=c;
                    else if (d=="tmRNAs") tmrnas[ctg]=c;
                    else if (d=="CRISPRs") crisprs[ctg]=c;
                }}
                END {{
                    for (ctg in cds) seen[ctg]=1;
                    for (ctg in trnas) seen[ctg]=1;
                    for (ctg in tmrnas) seen[ctg]=1;
                    for (ctg in crisprs) seen[ctg]=1;
                    for (ctg in seen)
                    print ctg, (ctg in cds?cds[ctg]:0), (ctg in trnas?trnas[ctg]:0), (ctg in tmrnas?tmrnas[ctg]:0), (ctg in crisprs?crisprs[ctg]:0);
                }}' {characterization_dir}/out_pharokka/{basename}/{basename}_cds_functions.tsv | sort -t $'\t' -k1,1;
            }} > {characterization_dir}/out_pharokka/{basename}/{basename}_cds_functions_wide.tsv""",
            rf"join -t $'\t' -1 1 -2 1 {characterization_dir}/out_pharokka/{basename}/{basename}_length_gc_cds_density_sorted.tsv {characterization_dir}/out_pharokka/{basename}/{basename}_cds_functions_wide.tsv | join -t $'\t' -1 1 -2 1 - {characterization_dir}/out_pharokka/{basename}/{basename}_top_hits_mash_inphared_sorted.tsv > {characterization_dir}/out_pharokka/{basename}/{basename}_final.tsv",
            rf"rm {characterization_dir}/out_pharokka/{basename}/{basename}_length_gc_cds_density_sorted.tsv {characterization_dir}/out_pharokka/{basename}/{basename}_top_hits_mash_inphared_sorted.tsv {characterization_dir}/out_pharokka/{basename}/{basename}_cds_functions_wide.tsv",
            #
            rf"""cut -f2- {characterization_dir}/out_taxmyphage/{basename}/Summary_taxonomy.tsv | awk 'BEGIN {{FS=OFS="\t"}} NR==1 {{$1="contig"; for (i=2; i<=NF; i++) $i="taxmyph_"$i}} 1' > {characterization_dir}/out_taxmyphage/{basename}/{basename}_taxmyphage_final.tsv""",
            #
            rf"""awk -F'\t' 'NR==1 {{print "contig\tphabox_length\tphabox_pred\tphabox_proportion\tphabox_pha_mer_score\tphabox_pha_mer_confidence\tphabox_lineage\tphabox_pha_gcn_score\tphabox_genus\tphabox_genus_cluster\tphabox_type\tphabox_pha_typ_score\tphabox_host\tphabox_cherry_score\tphabox_method\tphabox_host_ncbi_lineage\tphabox_host_gtdb_lineage\tphabox_phage_superkingdom\tphabox_phage_clade\tphabox_phage_kingdom\tphabox_phage_phylum\tphabox_phage_class\tphabox_phage_order\tphabox_phage_family\tphabox_phage_subfamily\tphabox_phage_genus\tphabox_host_domain_ncbi\tphabox_host_phylum_ncbi\tphabox_host_class_ncbi\tphabox_host_order_ncbi\tphabox_host_family_ncbi\tphabox_host_genus_ncbi\tphabox_host_species_ncbi"; next}} {{split($7,levels,";"); delete tax; for(i in levels){{split(levels[i],kv,":"); tax[kv[1]]=kv[2]}} split($16,levels2,";"); delete tax2; for(i in levels2){{split(levels2[i],kv2,"__"); tax2[kv2[1]]=kv2[2]}} print $0 "\t" tax["superkingdom"] "\t" tax["clade"] "\t" tax["kingdom"] "\t" tax["phylum"] "\t" tax["class"] "\t" tax["order"] "\t" tax["family"] "\t" tax["subfamily"] "\t" tax["genus"] "\t" tax2["d"] "\t" tax2["p"] "\t" tax2["c"] "\t" tax2["o"] "\t" tax2["f"] "\t" tax2["g"] "\t" tax2["s"]}}' {characterization_dir}/out_phabox/{basename}/final_prediction/final_prediction_summary.tsv > {characterization_dir}/out_phabox/{basename}/{basename}_phabox_final.tsv""",
            #
            rf"""awk 'BEGIN {{FS=OFS=","}} NR==1 {{$1="contig"; for (i=2; i<=NF; i++) $i="vhulk_"$i}} NR>1 {{$5=gensub(/_/, " ", "g", $5)}} 1' {characterization_dir}/out_vhulk/{basename}/{basename}_vHULK.csv | sed 's/,/\t/g' | sed 's/\r//' > {characterization_dir}/out_vhulk/{basename}/{basename}_vHULK_final.tsv""",
            #
            rf"""awk 'FNR==NR && /^>/ {{sub(/^>/, "", $1); sub(/ .*/, "", $1); contigs[$1]=0; next}} FNR>1 {{if ($2 in contigs) {{if (FILENAME ~ "card") card[$2]++; else if (FILENAME ~ "megares") megares[$2]++; else if (FILENAME ~ "ncbi") ncbi[$2]++; else if (FILENAME ~ "vfdb") vfdb[$2]++}}}} END {{print "contig\t", "abricate_card_hits\t", "abricate_megares_hits\t", "abricate_ncbi_hits\t", "abricate_vfdb_hits"; for (contig in contigs) printf "%s\t%d\t%d\t%d\t%d\n", contig, (card[contig] ? card[contig] : 0), (megares[contig] ? megares[contig] : 0), (ncbi[contig] ? ncbi[contig] : 0), (vfdb[contig] ? vfdb[contig] : 0)}}' {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes.fasta {characterization_dir}/out_abricate/{basename}/card_{basename}.abr {characterization_dir}/out_abricate/{basename}/megares_{basename}.abr {characterization_dir}/out_abricate/{basename}/ncbi_{basename}.abr {characterization_dir}/out_abricate/{basename}/vfdb_{basename}.abr > {characterization_dir}/out_abricate/{basename}/{basename}_abricate_final.tsv""",
            #
            rf"""join -t $'\t' <(sort -k1,1 {characterization_dir}/out_phabox/{basename}/{basename}_phabox_final.tsv) <(sort -k1,1 {characterization_dir}/out_pharokka/{basename}/{basename}_final.tsv) | join -t $'\t' - <(sort -k1,1 {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes_final.tsv) | join -t $'\t' - <(sort -k1,1 {characterization_dir}/out_abricate/{basename}/{basename}_abricate_final.tsv) | join -t $'\t' - <(sort -k1,1 {characterization_dir}/out_vhulk/{basename}/{basename}_vHULK_final.tsv) | join -t $'\t' - <(sort -k1,1 {characterization_dir}/out_taxmyphage/{basename}/{basename}_taxmyphage_final.tsv) | join -t $'\t' - <(sort -k1,1 {assembly_dir}/depth_coverage/{basename}_aln.tsv) | awk -v sample="{basename}" 'BEGIN {{OFS=FS="\t"}} NR==1 {{print "sample", $0}} NR>1 {{print sample, $0}}' > {characterization_dir}/intermediate/{basename}_characterization_tmp.tsv""",
            rf"""awk -F'\t' 'NR==1{{print $0 "\tpredicted_genus\tpredicted_species"; next}} {{
                delete count;
                delete species_by_genus;
                if($89!="")count[$89]++;
                if($60!="")count[$60]++;

                genus_from_14="";
                species92=$92;
                gsub(/^[[:space:]]+|[[:space:]]+$/, "", species92);
                species92_genus="";
                if(species92!=""){{
                    split(species92, sp92_parts, /[[:space:]]+/);
                    species92_genus=sp92_parts[1];
                }}
                if($14!="" && $14!="-"){{
                    n=split($14, entries, /\\n|\r/);
                    for(i=1; i<=n; i++){{
                        e=entries[i];
                        gsub(/^[[:space:]]+|[[:space:]]+$/, "", e);
                        if(e=="" || e=="-") continue;

                        if(e ~ /^genus:/){{
                            sub(/^genus:/, "", e);
                            split(e, parts, /[[:space:]]+/);
                            genus_from_14=parts[1];
                            break;
                        }}

                        if(e ~ /^species:/){{
                            species14=e;
                            sub(/^species:/, "", species14);
                            gsub(/^[[:space:]]+|[[:space:]]+$/, "", species14);
                            if(genus_from_14==""){{
                                split(species14, parts, /[[:space:]]+/);
                                genus_from_14=parts[1];
                            }}
                            split(species14, sp14_parts, /[[:space:]]+/);
                            sp14_genus=sp14_parts[1];
                            if(species_by_genus[sp14_genus]=="") species_by_genus[sp14_genus]=species14;
                        }}
                    }}
                }}
                if(genus_from_14!="")count[genus_from_14]++;

                predicted_genus="";
                max_count=0;
                for(val in count) if(count[val]>max_count){{predicted_genus=val; max_count=count[val]}}
                if(max_count<=1) predicted_genus="";

                predicted_species="";
                if(predicted_genus!=""){{
                    if(species92!="" && species92_genus==predicted_genus) predicted_species=species92;
                    else if(species_by_genus[predicted_genus]!="") predicted_species=species_by_genus[predicted_genus];
                }}
                delete species_by_genus;
                print $0 "\t" predicted_genus "\t" predicted_species
            }}' {characterization_dir}/intermediate/{basename}_characterization_tmp.tsv > {characterization_dir}/intermediate/{basename}_characterization_full.tsv""",
            rf"""rm {characterization_dir}/intermediate/{basename}_characterization_tmp.tsv"""
        ],
        7: [ # RESULTS
            f"echo \"Tool version (seqkit):\" $(seqkit version)",
            rf"mkdir -p {output_dir}/results/genomes/{basename}",
            rf"seqkit replace -p 'rotated.*' -r '' {characterization_dir}/out_pharokka/{basename}/{basename}_dnaapler_reoriented.fasta | seqkit seq -w 0 | seqkit split -i --by-id-prefix '{basename}_' -O {output_dir}/results/genomes/{basename}/",
            #
            rf"""awk -F'\t' '{{print $1","$2","$35","$36","$38","$39","$40","$41","$42","$110","$111","$106","$107","$109","$108","$112","$113","$75","$76","$77","$78","$79","$80","$81","$82","$83","$4","$12","$85","$86","$87","$88","$95","$96","$97","$98","$99","$100","$101","$102","$103","$60","$89","$14","$92","$114","$115}}' {characterization_dir}/intermediate/{basename}_characterization_full.tsv > {output_dir}/results/{basename}_final_results.csv"""
        ]
    }
    return commands_by_step.get(step_number, [])

def run_characterization(output_dir, input_files, task_params):
    characterization_dir = os.path.join(output_dir, "characterization_task")
    results_dir = os.path.join(output_dir, "results")

    # Ensure the characterization directory and subdirectories exist
    os.makedirs(characterization_dir, exist_ok=True)
    log_path = os.path.join(output_dir, "characterization.log")
    setup_logging(log_path)

    logger.info(f"Using folder: \"{characterization_dir}\"")

    for sub_dir in OUTPUT_CHAR_DIRS.values():
        os.makedirs(os.path.join(characterization_dir, sub_dir), exist_ok=True)

    for sub_dir in OUTPUT_RES_DIRS.values():
        os.makedirs(os.path.join(results_dir, sub_dir), exist_ok=True)

    # Start total timer
    total_start_time = time.time()

    # Define used environments per step
    envs = {
        1: {"env_name": "QPPL_env"},
        2: {"env_name": "QPPL_env2"},
        3: {"env_name": "QPPL_env3"},
        4: {"env_name": "vHULK"},
        5: {"env_name": "QPPL_env"},
        6: {"env_name": "QPPL_env"},
        7: {"env_name": "QPPL_env"}
    }

    # Process each input file
    for file in input_files:
        basename = re.sub(r'(\.fastq|\.fastq\.gz|\.fq|\.fq\.gz)$', '', file, flags=re.IGNORECASE)
        
        # Start task timer per file
        start_time = time.time()

        logger.info("=" * 60)
        logger.info(f"Starting characterization tasks for file: {file}")

        # Run commands in environments
        for step_number in sorted(envs.keys()):
            commands = generate_commands(characterization_dir, output_dir, file, task_params, basename, step_number)
            for command in commands:
                run_command(envs[step_number]['env_name'], command, file, step_number)

        # Log time taken for the file
        elapsed = time.time() - start_time
        logger.info(f"Finished characterization tasks for file: {file} (elapsed time: {elapsed:.2f} seconds)")
    
    # Calculate total elapsed time
    total_elapsed = time.time() - total_start_time
    logger.info(f"Total elapsed time for all tasks: {total_elapsed:.2f} seconds")