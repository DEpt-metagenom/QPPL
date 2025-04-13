import os
import subprocess
import re

def run_commands_in_env(env_name, commands, file, env_number):
    for command in commands:
        print(f"\n\033[1;32m[INFO] Running command in {env_name} on {file} (step {env_number}):\n  \033[0m$ {command}")
        result = subprocess.run(f"conda run -n {env_name} bash -c \"{command}\"", 
                                shell=True, text=True, capture_output=True)
        if result.returncode != 0:
            print(f"\033[1;31m[ERROR] Command failed: {command}\n\033[0m{result.stderr}")
        else:
            print(result.stdout)

def run_characterization(input_dir, output_dir, prefix, input_files, task_params):

    assembly_dir = os.path.join(output_dir, "assembly_task")

    characterization_dir = os.path.join(output_dir, "characterization_task")

    results_dir = os.path.join(output_dir, "results")

    sub_dirs_char = [
        "out_pharokka", "out_taxmyphage",
        "out_phabox", "out_vhulk",
        "out_abricate", "intermediate"
    ]

    sub_dirs_res = [
        "genomes"
    ]

    ### CREATE CHARACTERIZATION AND RESULTS DIRECTORY AND SUBDIRECTORIES
    if os.path.isdir(characterization_dir) and os.path.isdir(results_dir):
        print(f"\033[1;31m[INFO]\033[0m Detected existing folder: \"{characterization_dir}\"")
        print(f"\033[1;31m[INFO]\033[0m Detected existing folder: \"{results_dir}\"")
    else:
        print(f"\033[1;32m[INFO]\033[0m Created new folder: \"{characterization_dir}\"")
        print(f"\033[1;32m[INFO]\033[0m Created new folder: \"{results_dir}\"")
        os.makedirs(characterization_dir, exist_ok=True)
        os.makedirs(results_dir, exist_ok=True)

    for sub_dir in sub_dirs_char:
        path = os.path.join(characterization_dir, sub_dir)
        if not os.path.isdir(path):
            print(f"\033[1;32m[INFO]\033[0m Created new folder: \"{path}\"")
            os.makedirs(path, exist_ok=True)
        else:
            print(f"\033[1;31m[INFO]\033[0m Folder already exists: \"{path}\"")

    for sub_dir in sub_dirs_res:
        path = os.path.join(results_dir, sub_dir)
        if not os.path.isdir(path):
            print(f"\033[1;32m[INFO]\033[0m Created new folder: \"{path}\"")
            os.makedirs(path, exist_ok=True)
        else:
            print(f"\033[1;31m[INFO]\033[0m Folder already exists: \"{path}\"")

    ### LOOP OVER ALL INPUT FILES
    for file in input_files:
        basename = re.sub(r'(\.fastq|\.fastq\.gz|\.fq|\.fq\.gz)$', '', file)

        ### PRINT DEBUG INFORMATION
        print(f"\033[1;34m[DEBUG]\033[0m input_dir: {input_dir}")
        print(f"\033[1;34m[DEBUG]\033[0m output_dir: {output_dir}")
        print(f"\033[1;34m[DEBUG]\033[0m prefix: {prefix}")
        print(f"\033[1;34m[DEBUG]\033[0m file: {file}")
        print(f"\033[1;34m[DEBUG]\033[0m basename: {basename}")

        for key, value in task_params.get("characterization", {}).items():
            print(f"\033[1;34m[DEBUG]\033[0m {key}: {value}")

        ### RUN CHARACTERIZATION TOOLS
        envs = {
            1: {
                "env_name": "QPPL_env",
                "commands": [
                    rf"if [[ \$(grep -c \">\" {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes.fasta) -gt 1 ]]; then pharokka.py -i {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes.fasta -o {characterization_dir}/out_pharokka/{basename} -t {task_params['characterization']['pharokka_threads']} -m --dnaapler -p {basename} -d {task_params['characterization']['pharokka_db']} -g prodigal-gv; else pharokka.py -i {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes.fasta -o {characterization_dir}/out_pharokka/{basename} -t {task_params['characterization']['pharokka_threads']} --dnaapler -p {basename} -d {task_params['characterization']['pharokka_db']} -g prodigal-gv; fi",
                    rf"awk '/^LOCUS/ {{ if (out) close(out); out = \"{characterization_dir}/out_pharokka/{basename}/\" substr(\$2, 1) \".gbk\"; print \"Creating file:\", out; }} {{ print > out }}'  {characterization_dir}/out_pharokka/{basename}/{basename}.gbk",
                ]
            },
            2: {
                "env_name": "QPPL_env2",
                "commands": [
                    ##### INSTALL TAXMYPHAGE DB FIRST #####
                    rf"taxmyphage run -i {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes.fasta -t {task_params['characterization']['taxmyphage_threads']} -o {characterization_dir}/out_taxmyphage/{basename}",
                ]
            },
            3: {
                "env_name": "QPPL_env3",
                "commands": [
                    ##### INSTALL PHABOX DB FIRST #####
                    rf"phabox2 --task end_to_end --contigs {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes.fasta --outpth {characterization_dir}/out_phabox/{basename} --dbdir {task_params['characterization']['phabox_db']} --threads {task_params['characterization']['phabox_threads']}",
                ]
            },
            4: {
                "env_name": "vHULK",
                "commands": [
                    ##### INSTALL vHULK DB FIRST #####
                    rf"mkdir {task_params['characterization']['vhulk_location']}/{basename}_in",
                    rf"mkdir {task_params['characterization']['vhulk_location']}/{basename}",
                    rf"seqkit split -i {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes.fasta -O {task_params['characterization']['vhulk_location']}/{basename}_in --extension \".fasta\"",
                    rf"cd {task_params['characterization']['vhulk_location']} && python vHULK.py -i {basename}_in -o {basename} -t {task_params['characterization']['vhulk_threads']} && cd -",
                    rf"mv {task_params['characterization']['vhulk_location']}/{basename} {characterization_dir}/out_vhulk/{basename}",
                    rf"rm -r {task_params['characterization']['vhulk_location']}/{basename}_in",
                    rf"echo \"contig,pred_genus,score_genus,entropy_genus,pred_species,score_species,entropy_species\" > {characterization_dir}/out_vhulk/{basename}/{basename}_vHULK.csv",
                    rf"contigs=\$(ls {characterization_dir}/out_vhulk/{basename}/predictions/* | sed -n 's/.*part_\([^\.]*\)\.csv/\1/p') && for contig in \$contigs; do tail -n +2 {characterization_dir}/out_vhulk/{basename}/predictions/*\${{contig}}*.csv | awk -v sample=\"\${{contig}}\" '{{print sample \",\" \$0}}' >> {characterization_dir}/out_vhulk/{basename}/{basename}_vHULK.csv; done",
                ]
            },
            5: {
                "env_name": "QPPL_env",
                "commands": [
                    rf"mkdir {characterization_dir}/out_abricate/{basename}",
                    rf"abricate -db vfdb {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes.fasta --quiet > {characterization_dir}/out_abricate/{basename}/vfdb_{basename}.abr",
                    rf"abricate -db card {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes.fasta --quiet > {characterization_dir}/out_abricate/{basename}/card_{basename}.abr",
                    rf"abricate -db ncbi {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes.fasta --quiet > {characterization_dir}/out_abricate/{basename}/ncbi_{basename}.abr",
                    rf"abricate -db megares {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes.fasta --quiet > {characterization_dir}/out_abricate/{basename}/megares_{basename}.abr",
                    rf"abricate --summary {characterization_dir}/out_abricate/{basename}/vfdb_{basename}.abr > {characterization_dir}/out_abricate/{basename}/vfdb_{basename}.tab",
                    rf"abricate --summary {characterization_dir}/out_abricate/{basename}/card_{basename}.abr > {characterization_dir}/out_abricate/{basename}/card_{basename}.tab",
                    rf"abricate --summary {characterization_dir}/out_abricate/{basename}/ncbi_{basename}.abr > {characterization_dir}/out_abricate/{basename}/ncbi_{basename}.tab",
                    rf"abricate --summary {characterization_dir}/out_abricate/{basename}/megares_{basename}.abr > {characterization_dir}/out_abricate/{basename}/megares_{basename}.tab",
                ]
            },
            6: {
                "env_name": "QPPL_env",
                "commands": [
                    # CHECKV RESULTS
                    rf"awk 'BEGIN {{FS=OFS=\",\"}} NR==1 {{\$1=\"contig\"; \$2=\"vclust_cluster_id\"; \$NF=\"repres_score\"; for (i=3; i<=NF-1; i++) \$i=\"checkv_\"\$i}} 1' {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes.csv | sed 's/,/\t/g' > {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes_final.tsv",
                    # PHAROKKA RESULTS
                    rf"sort -t \$'\t' -k1,1 {characterization_dir}/out_pharokka/{basename}/{basename}_length_gc_cds_density.tsv > {characterization_dir}/out_pharokka/{basename}/{basename}_length_gc_cds_density_sorted.tsv",
                    rf"sort -t \$'\t' -k1,1 {characterization_dir}/out_pharokka/{basename}/{basename}_top_hits_mash_inphared.tsv > {characterization_dir}/out_pharokka/{basename}/{basename}_top_hits_mash_inphared_sorted.tsv",
                    rf"join -t \$'\t' -1 1 -2 1 {characterization_dir}/out_pharokka/{basename}/{basename}_length_gc_cds_density_sorted.tsv {characterization_dir}/out_pharokka/{basename}/{basename}_top_hits_mash_inphared_sorted.tsv > {characterization_dir}/out_pharokka/{basename}/{basename}_merged.tsv",
                    rf"awk 'BEGIN {{FS=OFS=\"\t\"}} NR==1 {{\$1=\"contig\"; for (i=2; i<=NF; i++) \$i=\"pharokka_\"\$i}} 1' {characterization_dir}/out_pharokka/{basename}/{basename}_merged.tsv > {characterization_dir}/out_pharokka/{basename}/{basename}_pharokka_final.tsv",
                    rf"rm {characterization_dir}/out_pharokka/{basename}/{basename}_length_gc_cds_density_sorted.tsv {characterization_dir}/out_pharokka/{basename}/{basename}_top_hits_mash_inphared_sorted.tsv {characterization_dir}/out_pharokka/{basename}/{basename}_merged.tsv",
                    # TAXMYPHAGE RESULTS
                    rf"cut -f2- {characterization_dir}/out_taxmyphage/{basename}/Summary_taxonomy.tsv | awk 'BEGIN {{FS=OFS=\"\t\"}} NR==1 {{\$1=\"contig\"; for (i=2; i<=NF; i++) \$i=\"taxmyph_\"\$i}} 1' > {characterization_dir}/out_taxmyphage/{basename}/{basename}_taxmyphage_final.tsv",
                    # PHABOX RESULTS
                    rf"awk -F'\t' 'NR==1 {{print \"contig\tphabox_length\tphabox_pred\tphabox_proportion\tphabox_pha_mer_score\tphabox_pha_mer_confidence\tphabox_lineage\tphabox_pha_gcn_score\tphabox_genus\tphabox_genus_cluster\tphabox_type\tphabox_pha_typ_score\tphabox_host\tphabox_cherry_score\tphabox_method\tphabox_host_ncbi_lineage\tphabox_host_gtdb_lineage\tphabox_phage_superkingdom\tphabox_phage_clade\tphabox_phage_kingdom\tphabox_phage_phylum\tphabox_phage_class\tphabox_phage_order\tphabox_phage_family\tphabox_phage_subfamily\tphabox_phage_genus\tphabox_host_domain_ncbi\tphabox_host_phylum_ncbi\tphabox_host_class_ncbi\tphabox_host_order_ncbi\tphabox_host_family_ncbi\tphabox_host_genus_ncbi\tphabox_host_species_ncbi\"; next}} {{split(\$7,levels,\";\"); delete tax; for(i in levels){{split(levels[i],kv,\":\"); tax[kv[1]]=kv[2]}} split(\$16,levels2,\";\"); delete tax2; for(i in levels2){{split(levels2[i],kv2,\"__\"); tax2[kv2[1]]=kv2[2]}} print \$0 \"\t\" tax[\"superkingdom\"] \"\t\" tax[\"clade\"] \"\t\" tax[\"kingdom\"] \"\t\" tax[\"phylum\"] \"\t\" tax[\"class\"] \"\t\" tax[\"order\"] \"\t\" tax[\"family\"] \"\t\" tax[\"subfamily\"] \"\t\" tax[\"genus\"] \"\t\" tax2[\"d\"] \"\t\" tax2[\"p\"] \"\t\" tax2[\"c\"] \"\t\" tax2[\"o\"] \"\t\" tax2[\"f\"] \"\t\" tax2[\"g\"] \"\t\" tax2[\"s\"]}}' {characterization_dir}/out_phabox/{basename}/final_prediction/final_prediction_summary.tsv > {characterization_dir}/out_phabox/{basename}/{basename}_phabox_final.tsv",
                    # VHULK RESULTS
                    rf"awk 'BEGIN {{FS=OFS=\",\"}} NR==1 {{\$1=\"contig\"; for (i=2; i<=NF; i++) \$i=\"vhulk_\"\$i}} NR>1 {{\$5=gensub(/_/, \" \", \"g\", \$5)}} 1' {characterization_dir}/out_vhulk/{basename}/{basename}_vHULK.csv | sed 's/,/\t/g' | sed 's/\r//' > {characterization_dir}/out_vhulk/{basename}/{basename}_vHULK_final.tsv",
                    # ABRICATE RESULTS
                    rf"awk 'FNR==NR && /^>/ {{sub(/^>/, \"\", \$1); sub(/ .*/, \"\", \$1); contigs[\$1]=0; next}} FNR>1 {{if (\$2 in contigs) {{if (FILENAME ~ \"card\") card[\$2]++; else if (FILENAME ~ \"megares\") megares[\$2]++; else if (FILENAME ~ \"ncbi\") ncbi[\$2]++; else if (FILENAME ~ \"vfdb\") vfdb[\$2]++}}}} END {{print \"contig\t\", \"abricate_card_hits\t\", \"abricate_megares_hits\t\", \"abricate_ncbi_hits\t\", \"abricate_vfdb_hits\"; for (contig in contigs) printf \"%s\t%d\t%d\t%d\t%d\n\", contig, (card[contig] ? card[contig] : 0), (megares[contig] ? megares[contig] : 0), (ncbi[contig] ? ncbi[contig] : 0), (vfdb[contig] ? vfdb[contig] : 0)}}' {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes.fasta {characterization_dir}/out_abricate/{basename}/card_{basename}.abr {characterization_dir}/out_abricate/{basename}/megares_{basename}.abr {characterization_dir}/out_abricate/{basename}/ncbi_{basename}.abr {characterization_dir}/out_abricate/{basename}/vfdb_{basename}.abr > {characterization_dir}/out_abricate/{basename}/{basename}_abricate_final.tsv",
                    # CONCATENATED RESULTS
                    rf"join -t \$'\t' <(sort -k1,1 {characterization_dir}/out_phabox/{basename}/{basename}_phabox_final.tsv) <(sort -k1,1 {characterization_dir}/out_pharokka/{basename}/{basename}_pharokka_final.tsv) | join -t \$'\t' - <(sort -k1,1 {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes_final.tsv) | join -t \$'\t' - <(sort -k1,1 {characterization_dir}/out_abricate/{basename}/{basename}_abricate_final.tsv) | join -t \$'\t' - <(sort -k1,1 {characterization_dir}/out_vhulk/{basename}/{basename}_vHULK_final.tsv) | join -t \$'\t' - <(sort -k1,1 {characterization_dir}/out_taxmyphage/{basename}/{basename}_taxmyphage_final.tsv) | awk -v sample=\"{basename}\" 'BEGIN {{OFS=FS=\"\t\"}} NR==1 {{print \"sample\", \$0}} NR>1 {{print sample, \$0}}' > {characterization_dir}/intermediate/{basename}_characterization_tmp.tsv",
                    rf"awk -F'\t' 'NR==1{{print \$0 \"\tpredicted_genus\tpredicted_species\"; next}} {{delete count; if(\$85!=\"\")count[\$85]++; if(\$56!=\"\")count[\$56]++; if(\$33!=\"\")count[\$33]++; predicted_genus=\"\"; max_count=0; for(val in count) if(count[val]>max_count){{predicted_genus=val; max_count=count[val]}} if(max_count<=1) predicted_genus=\"\"; predicted_species=\"\"; if(\$88!=\"\"&&\$34!=\"\"&&\$88==\$34) predicted_species=\$88; print \$0 \"\t\" predicted_genus \"\t\" predicted_species}}' {characterization_dir}/intermediate/{basename}_characterization_tmp.tsv > {characterization_dir}/intermediate/{basename}_characterization_full.tsv",
                    rf"rm {characterization_dir}/intermediate/{basename}_characterization_tmp.tsv",
                ]
            },
            7: {
                "env_name": "QPPL_env",
                "commands": [
                    # MOVE GENOMES TO RESULTS DIRECTORY
                    rf"mkdir -p {results_dir}/genomes/{basename}",
                    rf"seqkit replace -p 'rotated.*' -r '' {characterization_dir}/out_pharokka/{basename}/{basename}_dnaapler_reoriented.fasta | seqkit seq -w 0 | seqkit split -i --by-id-prefix '{basename}_' -O {results_dir}/genomes/{basename}/",
                    # CREATE FINAL TABLE
                    rf"awk -F'\t' '{{print \$1\",\"\$2\",\"\$35\",\"\$36\",\"\$50\",\"\$38\",\"\$51\",\"\$52\",\"\$53\",\"\$55\",\"\$4\",\"\$12\",\"\$81\",\"\$82\",\"\$83\",\"\$84\",\"\$91\",\"\$92\",\"\$93\",\"\$94\",\"\$95\",\"\$96\",\"\$97\",\"\$98\",\"\$99\",\"\$102\",\"\$103}}' {characterization_dir}/intermediate/{basename}_characterization_full.tsv > {results_dir}/{basename}_final_results.csv"
                ]
            }
        }

        for step_number in sorted(envs.keys()):
            env = envs[step_number]
            run_commands_in_env(env["env_name"], env["commands"], file, step_number)