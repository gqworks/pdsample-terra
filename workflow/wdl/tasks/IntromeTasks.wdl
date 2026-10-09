version 1.0

# Native WDL port of Introme v2 infer mode.
# Pin: https://github.com/CCICB/introme/commit/d6da48bc84abda66f588470f79abcf1244acfb21
# (branch nextflow_fns as of 2026-08-27). Squirls is not called: main.nf comments it out.
# ensemble_train is not ported.
#
# Process -> image (CPU profile; GPU images are selected by the workflow):
#   data_preprocessing, quality_filter  gabyou/data_prepocessing:v2.0
#   variant_info                        gabyou/variant_info:v3.2
#   spliceai                            headoncollusion/spliceai:v1.4-cpu
#   mmsplice                            headoncollusion/mmsplice:v2.2-cpu
#   pangolin                            headoncollusion/pangolin:v1.3-cpu
#   spip                                headoncollusion/spip:v1.1
#   spliceogen                          headoncollusion/spliceogen:v2.1
#   introme_functions, splicing_anno,
#   ensemble_infer                      headoncollusion/ag_check:v1.6

task FetchIntromeAssets {
  input {
    File? introme_assets_tar
    String introme_commit = "d6da48bc84abda66f588470f79abcf1244acfb21"
    String archive_url = "https://github.com/CCICB/introme/archive/d6da48bc84abda66f588470f79abcf1244acfb21.tar.gz"
    String genome_build = "hg38"
    String deprecated_introme_docker = ""
    String force_rerun_token = ""
    String python_docker = "python:3.11.11-slim-bookworm"
    Int preemptible = 2
  }
  String supplied_tar = if defined(introme_assets_tar) then select_first([introme_assets_tar]) else ""
  command <<<
    set -euo pipefail
    echo "force_rerun_token=~{force_rerun_token}" >&2
    echo "introme_commit=~{introme_commit}" >&2
    if [[ -n "~{deprecated_introme_docker}" ]]; then
      echo "deprecated introme_docker is ignored: ~{deprecated_introme_docker}" >&2
    fi
    if [[ "~{genome_build}" != "hg38" ]]; then
      echo "Introme v2 WDL port stages hg38 assets only. Upstream nextflow_fns marks hg19 untested." >&2
      exit 1
    fi
    pin="~{introme_commit}"
    mkdir -p src/unpack stage/tomls stage/branchpointer stage/regions stage/U12 stage/AG_check stage/ESE stage/ml stage/models
    if [[ -n "~{supplied_tar}" ]]; then
      cp "~{supplied_tar}" src/in.tar.gz
    else
      python3 -c "import urllib.request; urllib.request.urlretrieve('~{archive_url}', 'src/in.tar.gz')"
    fi
    tar -xzf src/in.tar.gz -C src/unpack
    if [[ -f src/unpack/COMMIT ]]; then
      got="$(tr -d '[:space:]' < src/unpack/COMMIT)"
      if [[ "$got" != "$pin" ]]; then
        echo "assets COMMIT ${got} does not match pin ${pin}" >&2
        exit 1
      fi
      tar -C src/unpack -czf introme-assets.tar .
    else
      shopt -s nullglob
      tops=(src/unpack/introme-*)
      if [[ ${#tops[@]} -ne 1 || ! -d "${tops[0]}" ]]; then
        echo "Expected a GitHub archive root introme-<sha>/ or a slim tar with COMMIT" >&2
        exit 1
      fi
      root="${tops[0]}"
      sha="$(basename "$root")"
      sha="${sha#introme-}"
      if [[ "$sha" != "$pin" ]]; then
        echo "archive directory $(basename "$root") does not match pin ${pin}" >&2
        exit 1
      fi
      copy() {
        if [[ ! -f "$1" ]]; then
          echo "Pinned archive is missing $1" >&2
          exit 1
        fi
        cp -L "$1" "$2"
      }
      copy "$root/nextflow/assets/chrRename.tsv" stage/chrRename.tsv
      copy "$root/nextflow/assets/conf_pre_anno.lua" stage/conf_pre_anno.lua
      copy "$root/nextflow/assets/conf_ensemble.lua" stage/conf_ensemble.lua
      copy "$root/nextflow/assets/introme_annotate.vcf" stage/introme_annotate.vcf
      copy "$root/nextflow/assets/tomls/gencode.hg38.toml" stage/tomls/gencode.hg38.toml
      copy "$root/nextflow/assets/tomls/annotate.hg38.toml" stage/tomls/annotate.hg38.toml
      copy "$root/nextflow/assets/tomls/vcfanno_splicing_ensemble.toml" stage/tomls/vcfanno_splicing_ensemble.toml
      copy "$root/nextflow/assets/tomls/introme.toml" stage/tomls/introme.toml
      copy "$root/nextflow/assets/branchpointer/branchpointer.hg38.bed.gz" stage/branchpointer/branchpointer.hg38.bed.gz
      copy "$root/nextflow/assets/branchpointer/branchpointer.hg38.bed.gz.tbi" stage/branchpointer/branchpointer.hg38.bed.gz.tbi
      copy "$root/nextflow/assets/regions/regions.hg38.bed.gz" stage/regions/regions.hg38.bed.gz
      copy "$root/nextflow/assets/regions/regions.hg38.bed.gz.tbi" stage/regions/regions.hg38.bed.gz.tbi
      copy "$root/nextflow/assets/U12/U12.hg38.bed.gz" stage/U12/U12.hg38.bed.gz
      copy "$root/nextflow/assets/U12/U12.hg38.bed.gz.tbi" stage/U12/U12.hg38.bed.gz.tbi
      copy "$root/AG_check/AG_check.py" stage/AG_check/AG_check.py
      copy "$root/ESE/scoring.py" stage/ESE/scoring.py
      copy "$root/ESE/motifs.py" stage/ESE/motifs.py
      copy "$root/ESE/variants.py" stage/ESE/variants.py
      copy "$root/ESE/ESEfinder_motif_source.py" stage/ESE/ESEfinder_motif_source.py
      copy "$root/ESE/RCRUNCH_motif_source.py" stage/ESE/RCRUNCH_motif_source.py
      copy "$root/nextflow/modules/ML/main.py" stage/ml/main.py
      copy "$root/nextflow/modules/ML/infer.py" stage/ml/infer.py
      copy "$root/nextflow/modules/ML/train.py" stage/ml/train.py
      copy "$root/nextflow/modules/ML/utils.py" stage/ml/utils.py
      copy "$root/nextflow/modules/ML/pipeline_constants.py" stage/ml/pipeline_constants.py
      copy "$root/nextflow/assets/models/20042026/RandomForest_all_info_introme_train_run_model.pkl" stage/models/RandomForest_all_info_introme_train_run_model.pkl
      copy "$root/nextflow/assets/models/20042026/RandomForest_all_info_introme_train_run_columns.json" stage/models/RandomForest_all_info_introme_train_run_columns.json
      printf '%s\n' "$pin" > stage/COMMIT
      tar -C stage -czf introme-assets.tar .
    fi
    printf '%s\n' "$pin" > commit.txt
  >>>
  output {
    File assets_tar = "introme-assets.tar"
    String commit_sha = read_string("commit.txt")
  }
  runtime {
    docker: python_docker
    preemptible: preemptible
    memory: "2 GiB"
    disks: "local-disk 20 HDD"
    cpu: 1
  }
}

task DataPreprocessing {
  input {
    File vcf
    File vcf_index
    File ref_fasta
    File ref_fasta_index
    File gtf
    File assets_tar
    File? bed
    Boolean use_bed = true
    String prefix = "cohort"
    String force_rerun_token = ""
    String docker = "gabyou/data_prepocessing:v2.0"
    Int preemptible = 2
    Int cpu = 4
    Int memory_gb = 8
  }
  String bed_path = if defined(bed) then select_first([bed]) else ""
  String gtf_prefix = if basename(gtf) != basename(gtf, ".gtf.gz") then basename(gtf, ".gtf.gz") else basename(gtf, ".gtf")
  Int disk_gb = ceil(size(vcf, "GiB") + size(ref_fasta, "GiB") + size(gtf, "GiB") + size(assets_tar, "GiB") + 20)
  command <<<
    set -euo pipefail
    echo "force_rerun_token=~{force_rerun_token}" >&2
    tar -xzf ~{assets_tar}
    echo "introme_commit=$(cat COMMIT)" >&2
    mkdir -p ref
    ln -s "$(readlink -f ~{ref_fasta})" ref/ref.fa
    ln -s "$(readlink -f ~{ref_fasta_index})" ref/ref.fa.fai
    # Faithful port of nextflow/modules/data_preprocessing.nf.
    # publishDir and docker runOptions -u are Cromwell's job, not part of the script.
    CHROMS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
    gtf_path="~{gtf_prefix}.sorted.gtf.gz"
    input_gtf="~{gtf}"
    input_vcf="~{vcf}"
    echo $(date +%x_%r) $(bcftools view -H "$input_vcf" | wc -l) 'variants prior to subsetting'
    bcftools annotate --rename-chrs chrRename.tsv -Ou "$input_vcf" \
        | bcftools sort -Ou \
        | bcftools view -t "$CHROMS" -Ou \
        | bcftools norm -f ref/ref.fa -c e -m-both -Ou \
        | bcftools norm -d exact -Ou \
        | bcftools sort -Oz -o ~{prefix}.sorted.norm.vcf.gz
    if tabix -f -p gff "$input_gtf"; then
        cp -L "$input_gtf" "$gtf_path"
    else
        echo $(date +%x_%r) 'GTF file not in BGZF format - sorting and bgzipping for tabix indexing'
        (
            gunzip -c "$input_gtf" | grep '^#' || true
            gunzip -c "$input_gtf" | grep -v '^#' |
                sort -k1,1 -k4,4n -k5,5n -s
        ) | bgzip > "$gtf_path"
        tabix -f "$gtf_path"
    fi
    rm -f "$gtf_path.tbi"
    bed_path="~{bed_path}"
    if [ "~{use_bed}" != "true" ] || [ -z "$bed_path" ]; then
        echo $(date +%x_%r) 'No BED file provided - Beginning subsetting to GTF regions'
        bedtools intersect -header -u -a ~{prefix}.sorted.norm.vcf.gz -b "$gtf_path" | bgzip > ~{prefix}.subset.vcf.gz
    else
        echo $(date +%x_%r) 'BED file provided - Beginning subsetting to genomic regions of interest'
        bedtools intersect -header -u -a ~{prefix}.sorted.norm.vcf.gz -b "$bed_path" | bgzip > ~{prefix}.subset.vcf.gz
    fi
    tabix -p vcf ~{prefix}.subset.vcf.gz
    variant_count=$(bcftools view -H ~{prefix}.subset.vcf.gz | wc -l | tr -d ' ')
    if [ "$variant_count" -eq "0" ]; then
        if [ -z "$bed_path" ]; then
            echo $(date +%x_%r) 'No variants were in regions of interest - is your GTF file restricted to certain regions?'
        else
            echo $(date +%x_%r) 'No variants were in regions of interest - perhaps expand your BED file to more regions of interest'
        fi
        exit 1
    else
        # Upstream message says "prior" but prints the post-subset count.
        echo $(date +%x_%r) "$variant_count variants prior to subsetting"
    fi
  >>>
  output {
    File subset_vcf = "~{prefix}.subset.vcf.gz"
    File subset_vcf_index = "~{prefix}.subset.vcf.gz.tbi"
    File sorted_gtf = "~{gtf_prefix}.sorted.gtf.gz"
  }
  runtime {
    docker: docker
    preemptible: preemptible
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: cpu
    bootDiskSizeGb: 20
  }
}

task QualityFilter {
  input {
    File vcf
    File vcf_index
    String prefix = "cohort"
    Int min_qual = 200
    Int min_dp = 20
    Int min_ad = 5
    String force_rerun_token = ""
    String docker = "gabyou/data_prepocessing:v2.0"
    Int preemptible = 2
    Int cpu = 2
    Int memory_gb = 4
  }
  Int disk_gb = ceil(size(vcf, "GiB") + 20)
  command <<<
    set -euo pipefail
    echo "force_rerun_token=~{force_rerun_token}" >&2
    # Faithful port of nextflow/modules/quality_filter.nf. Default off in params.json.
    input_vcf="~{vcf}"
    has_format_dp=$(bcftools view -h ${input_vcf} | grep -q '^##FORMAT=<ID=DP,' && echo true || echo false)
    has_format_ad=$(bcftools view -h ${input_vcf} | grep -q '^##FORMAT=<ID=AD,' && echo true || echo false)
    expr="(FILTER='PASS' || FILTER='.') && (QUAL>=~{min_qual} || QUAL='.')"
    if [[ "$has_format_dp" == true ]]; then
        expr="$expr && MAX(FORMAT/DP[*])>=~{min_dp}"
    else
        echo "WARNING: FORMAT/DP not found in VCF header; skipping DP filter" >&2
    fi
    if [[ "$has_format_ad" == true ]]; then
        expr="$expr && MAX(FORMAT/AD[*:1])>=~{min_ad}"
    else
        echo "WARNING: FORMAT/AD not found in VCF header; skipping AD filter" >&2
    fi
    echo "Using bcftools filter expression: $expr" >&2
    bcftools filter --threads $(getconf _NPROCESSORS_ONLN) \
        -i "$expr" \
        ${input_vcf} \
        -Oz -o ~{prefix}.quality_filter.vcf.gz
    tabix -p vcf ~{prefix}.quality_filter.vcf.gz
    variant_count=$(bcftools view -H ~{prefix}.quality_filter.vcf.gz | wc -l | tr -d ' ')
    echo "$variant_count" > ~{prefix}.quality_filter.variant_count.txt
    echo $(date +%x_%r) 'Quality filtering complete -' $variant_count 'variants remaining'
    if [[ $variant_count == 0 ]]; then
        echo $(date +%x_%r) 'No variants passed quality filtering - perhaps rerun with -q (triggers no quality filtering)'
        exit 1
    fi
  >>>
  output {
    File quality_vcf = "~{prefix}.quality_filter.vcf.gz"
    File quality_vcf_index = "~{prefix}.quality_filter.vcf.gz.tbi"
    File variant_count = "~{prefix}.quality_filter.variant_count.txt"
  }
  runtime {
    docker: docker
    preemptible: preemptible
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: cpu
    bootDiskSizeGb: 20
  }
}

task VariantInfo {
  input {
    File vcf
    File vcf_index
    File gtf
    File assets_tar
    String prefix = "cohort"
    String genome_build = "hg38"
    Boolean filter_rare = true
    Float allele_frequency = 0.01
    String af_field = ""
    String force_rerun_token = ""
    String docker = "gabyou/variant_info:v3.2"
    Int preemptible = 2
    Int cpu = 4
    Int memory_gb = 8
  }
  Int disk_gb = ceil(size(vcf, "GiB") + size(gtf, "GiB") + size(assets_tar, "GiB") + 20)
  command <<<
    set -euo pipefail
    echo "force_rerun_token=~{force_rerun_token}" >&2
    tar -xzf ~{assets_tar}
    echo "introme_commit=$(cat COMMIT)" >&2
    # Faithful port of nextflow/modules/variant_info.nf (annotation + rmanno + stripped).
    # The gnomAD_PM_AF filter in that module is commented out. The block below is the
    # cost-control switch (introme_filter_rare). It runs only when an AF INFO field
    # is already on the VCF, because gencode.*.toml leaves gnomAD annotation commented out.
    toml="tomls/gencode.~{genome_build}.toml"
    if [[ ! -f "$toml" ]]; then
      echo "No $toml in the pinned assets. This port stages hg38 only." >&2
      exit 1
    fi
    # Copy the GTF into the writable work dir. Cromwell inputs are often read-only,
    # and tabix writes the index beside the file.
    cp -L "~{gtf}" sorted.annotation.gtf.gz
    tabix -f sorted.annotation.gtf.gz
    sed "s|REPLACE_gencode_file|${PWD}/sorted.annotation.gtf.gz|" "$toml" > temp.toml
    vcfanno -p $(getconf _NPROCESSORS_ONLN) -lua conf_pre_anno.lua temp.toml ~{vcf} | bgzip > ~{prefix}.variant_info.vcf.gz
    tabix -p vcf ~{prefix}.variant_info.vcf.gz
    ~{true="FILTER_RARE=1" false="FILTER_RARE=0" filter_rare}
    if [[ "$FILTER_RARE" == "1" ]]; then
      header="$(bcftools view -h ~{prefix}.variant_info.vcf.gz)"
      af_field="~{af_field}"
      if [[ -z "$af_field" ]]; then
        for candidate in gnomAD_PM_AF gnomAD_AF AF_popmax gnomADg_AF_popmax gnomADe_AF_popmax gnomad_AF_popmax; do
          if grep -q "^##INFO=<ID=${candidate}," <<<"$header"; then
            af_field="$candidate"
            break
          fi
        done
      fi
      if [[ -n "$af_field" ]]; then
        if ! grep -q "^##INFO=<ID=${af_field}," <<<"$header"; then
          echo "AF field ${af_field} is not in the VCF header" >&2
          exit 1
        fi
        echo "Rare filter ${af_field} <= ~{allele_frequency} (missing AF kept)" >&2
        bcftools filter --threads $(getconf _NPROCESSORS_ONLN) \
          -i "(${af_field}<=~{allele_frequency} || ${af_field}='.')" \
          -Oz -o ~{prefix}.variant_info.rare.vcf.gz ~{prefix}.variant_info.vcf.gz
        mv ~{prefix}.variant_info.rare.vcf.gz ~{prefix}.variant_info.vcf.gz
        tabix -f -p vcf ~{prefix}.variant_info.vcf.gz
        rare_count="$(bcftools view -H ~{prefix}.variant_info.vcf.gz | wc -l | tr -d ' ')"
        echo "$rare_count variants after allele-frequency filter" >&2
        if [[ "$rare_count" == "0" ]]; then
          echo "No variants passed allele_frequency ~{allele_frequency}. Raise it or set introme_filter_rare false." >&2
          exit 1
        fi
      else
        echo "WARNING: introme_filter_rare is set but no gnomAD INFO field is on the VCF. Upstream v2 comments out gnomAD annotation in gencode toml, so the filter is a no-op unless the input already carries gnomAD_PM_AF or a related INFO field. Cohort INFO/AF is not used." >&2
      fi
    fi
    gunzip -k ~{prefix}.variant_info.vcf.gz
    bcftools view -h ~{prefix}.variant_info.vcf.gz > ~{prefix}.variant_info.filtered.rmanno.vcf
    grep -v "^#" ~{prefix}.variant_info.vcf | awk '{$8="."; print }' OFS='\t' >> ~{prefix}.variant_info.filtered.rmanno.vcf
    bcftools view -G ~{prefix}.variant_info.vcf > ~{prefix}.variant_info.vcf.stripped
  >>>
  output {
    File variant_info_vcf = "~{prefix}.variant_info.vcf.gz"
    File variant_info_vcf_index = "~{prefix}.variant_info.vcf.gz.tbi"
    File variant_info_stripped = "~{prefix}.variant_info.vcf.stripped"
    File variant_info_rmanno = "~{prefix}.variant_info.filtered.rmanno.vcf"
  }
  runtime {
    docker: docker
    preemptible: preemptible
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: cpu
    bootDiskSizeGb: 20
  }
}

task SpliceAI {
  input {
    File vcf
    File ref_fasta
    File ref_fasta_index
    String prefix = "cohort"
    String genome_build = "hg38"
    String spliceai_db = "gencode.v44.annotation.txt.gz"
    File? spliceai_db_file
    Int distance = 1000
    Int mask = 0
    Boolean use_gpu = false
    String gpu_type = "nvidia-tesla-t4"
    Int gpu_count = 0
    String nvidia_driver_version = "535.183.01"
    String force_rerun_token = ""
    String docker = "headoncollusion/spliceai:v1.4-cpu"
    Int preemptible = 1
    Int cpu = 8
    Int memory_gb = 8
  }
  String db_path = if defined(spliceai_db_file) then select_first([spliceai_db_file]) else ""
  Int disk_gb = ceil(size(vcf, "GiB") + size(ref_fasta, "GiB") + 30)
  Int boot_gb = if use_gpu then 50 else 30
  command <<<
    set -euo pipefail
    echo "force_rerun_token=~{force_rerun_token}" >&2
    mkdir -p ref
    ln -s "$(readlink -f ~{ref_fasta})" ref/ref.fa
    ln -s "$(readlink -f ~{ref_fasta_index})" ref/ref.fa.fai
    # Faithful port of nextflow/modules/spliceai/spliceai.nf.
    # Optional spliceai_db_file skips the wget (stage it to avoid the runtime download).
    export TF_FORCE_GPU_ALLOW_GROWTH=true
    db_path="~{db_path}"
    if [[ -n "$db_path" ]]; then
      db="$db_path"
    else
      wget "https://compbio.ccia.org.au/introme/files/~{genome_build}/~{spliceai_db}" --no-check-certificate -O "~{spliceai_db}"
      db="~{spliceai_db}"
    fi
    touch ~{prefix}.spliceai.vcf
    ls -lah
    spliceai -I ~{vcf} -O ~{prefix}.spliceai.vcf -R ref/ref.fa -A "$db" -D ~{distance} -M ~{mask} 1>./log
    sed -i -E 's/([0-9]+),(NM_|ENST|ENSG)/\1\&\2/g' ~{prefix}.spliceai.vcf
  >>>
  output {
    File spliceai_vcf = "~{prefix}.spliceai.vcf"
  }
  runtime {
    docker: docker
    preemptible: preemptible
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: cpu
    bootDiskSizeGb: boot_gb
    gpuCount: gpu_count
    gpuType: gpu_type
    nvidiaDriverVersion: nvidia_driver_version
  }
}

task MMSplice {
  input {
    File vcf
    File ref_fasta
    File ref_fasta_index
    File gtf
    String prefix = "cohort"
    Boolean use_gpu = false
    String gpu_type = "nvidia-tesla-t4"
    Int gpu_count = 0
    String nvidia_driver_version = "535.183.01"
    String force_rerun_token = ""
    String docker = "headoncollusion/mmsplice:v2.2-cpu"
    Int preemptible = 1
    Int cpu = 8
    Int memory_gb = 8
  }
  Int disk_gb = ceil(size(vcf, "GiB") + size(ref_fasta, "GiB") + size(gtf, "GiB") + 20)
  Int boot_gb = if use_gpu then 50 else 30
  command <<<
    set -euo pipefail
    echo "force_rerun_token=~{force_rerun_token}" >&2
    mkdir -p ref
    ln -s "$(readlink -f ~{ref_fasta})" ref/ref.fa
    ln -s "$(readlink -f ~{ref_fasta_index})" ref/ref.fa.fai
    # Faithful port of nextflow/modules/mmsplice/mmsplice.nf.
    # run_mmsplice.py is on PATH inside the image. The repo copy is not mounted.
    pwd
    ls
    export TF_FORCE_GPU_ALLOW_GROWTH=true
    run_mmsplice.py --vcf ~{vcf} --fasta ref/ref.fa --gtf ~{gtf} --output ~{prefix}.mmsplice.vcf
  >>>
  output {
    File mmsplice_vcf = "~{prefix}.mmsplice.vcf"
  }
  runtime {
    docker: docker
    preemptible: preemptible
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: cpu
    bootDiskSizeGb: boot_gb
    gpuCount: gpu_count
    gpuType: gpu_type
    nvidiaDriverVersion: nvidia_driver_version
  }
}

task Pangolin {
  input {
    File vcf
    File ref_fasta
    File ref_fasta_index
    String prefix = "cohort"
    String genome_build = "hg38"
    String pangolin_db = "gencode.v44.annotation.db"
    File? pangolin_db_file
    Boolean use_gpu = false
    String gpu_type = "nvidia-tesla-t4"
    Int gpu_count = 0
    String nvidia_driver_version = "535.183.01"
    String force_rerun_token = ""
    String docker = "headoncollusion/pangolin:v1.3-cpu"
    Int preemptible = 1
    Int cpu = 8
    Int memory_gb = 8
  }
  String db_path = if defined(pangolin_db_file) then select_first([pangolin_db_file]) else ""
  Int disk_gb = ceil(size(vcf, "GiB") + size(ref_fasta, "GiB") + 30)
  Int boot_gb = if use_gpu then 50 else 30
  command <<<
    set -euo pipefail
    echo "force_rerun_token=~{force_rerun_token}" >&2
    mkdir -p ref
    ln -s "$(readlink -f ~{ref_fasta})" ref/ref.fa
    ln -s "$(readlink -f ~{ref_fasta_index})" ref/ref.fa.fai
    # Faithful port of nextflow/modules/pangolin/pangolin.nf.
    # The annotation DB is about 300 MB. Pass pangolin_db_file to skip wget.
    db_path="~{db_path}"
    if [[ -n "$db_path" ]]; then
      db="$db_path"
    else
      wget "https://compbio.ccia.org.au/introme/files/~{genome_build}/~{pangolin_db}" --no-check-certificate -O "~{pangolin_db}"
      db="~{pangolin_db}"
    fi
    pangolin ~{vcf} ref/ref.fa "$db" ~{prefix}.pangolin
  >>>
  output {
    File pangolin_vcf = "~{prefix}.pangolin.vcf"
  }
  runtime {
    docker: docker
    preemptible: preemptible
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: cpu
    bootDiskSizeGb: boot_gb
    gpuCount: gpu_count
    gpuType: gpu_type
    nvidiaDriverVersion: nvidia_driver_version
  }
}

task Spip {
  input {
    File vcf
    String prefix = "cohort"
    String genome_build = "hg38"
    String force_rerun_token = ""
    String docker = "headoncollusion/spip:v1.1"
    Int preemptible = 2
    Int cpu = 1
    Int memory_gb = 4
  }
  Int disk_gb = ceil(size(vcf, "GiB") + 20)
  command <<<
    set -euo pipefail
    echo "force_rerun_token=~{force_rerun_token}" >&2
    # Faithful port of nextflow/modules/spip/spip.nf.
    # The historical containerOptions bind-mount of outdir onto /data is commented out
    # upstream and is not used. Cromwell localizes the VCF and this command reads that path.
    # docker runOptions -u is not set; Cromwell collects outputs as the container user.
    ls
    pwd
    Rscript /SPiP/SPiPv2.1_main.r -I ~{vcf} -O ~{prefix}.spip.vcf -g ~{genome_build} --VCF
  >>>
  output {
    File spip_vcf = "~{prefix}.spip.vcf"
  }
  runtime {
    docker: docker
    preemptible: preemptible
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: cpu
    bootDiskSizeGb: 20
  }
}

task Spliceogen {
  input {
    File vcf
    File ref_fasta
    File ref_fasta_index
    File gtf
    String prefix = "cohort"
    String force_rerun_token = ""
    String docker = "headoncollusion/spliceogen:v2.1"
    Int preemptible = 2
    Int cpu = 3
    Int memory_gb = 4
  }
  Int disk_gb = ceil(size(vcf, "GiB") + size(ref_fasta, "GiB") + size(gtf, "GiB") + 20)
  command <<<
    set -euo pipefail
    echo "force_rerun_token=~{force_rerun_token}" >&2
    mkdir -p ref
    ln -s "$(readlink -f ~{ref_fasta})" ref/ref.fa
    ln -s "$(readlink -f ~{ref_fasta_index})" ref/ref.fa.fai
    # Faithful port of nextflow/modules/spliceogen/spliceogen.nf.
    # Nextflow sets containerOptions '--entrypoint ""'. Cromwell launches with bash,
    # which replaces the image entrypoint the same way.
    # Upstream memory is 1 GB at 3 CPUs. That is below the GCP custom-machine
    # minimum (~0.9 GiB/vCPU), so this task asks for 4 GiB.
    echo pwd
    pwd
    echo ls
    ls
    ORIG_DIR="$(pwd)"
    cd /Spliceogen
    ./RUN.sh -input "$(realpath "$ORIG_DIR/~{vcf}")" \
             -fasta "$(realpath "$ORIG_DIR/ref/ref.fa")" \
             -gtf "$(realpath "$ORIG_DIR/~{gtf}")"
    ls output
    mv output/*_out.txt "$ORIG_DIR/~{prefix}.spliceogen.tsv"
  >>>
  output {
    File spliceogen_tsv = "~{prefix}.spliceogen.tsv"
  }
  runtime {
    docker: docker
    preemptible: preemptible
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: cpu
    bootDiskSizeGb: 20
  }
}

task IntromeFunctions {
  input {
    File variant_info_stripped
    File ref_fasta
    File ref_fasta_index
    File assets_tar
    String prefix = "cohort"
    String force_rerun_token = ""
    String docker = "headoncollusion/ag_check:v1.6"
    Int preemptible = 2
    Int cpu = 2
    Int memory_gb = 8
  }
  Int disk_gb = ceil(size(variant_info_stripped, "GiB") + size(ref_fasta, "GiB") + size(assets_tar, "GiB") + 20)
  command <<<
    set -euo pipefail
    echo "force_rerun_token=~{force_rerun_token}" >&2
    tar -xzf ~{assets_tar}
    echo "introme_commit=$(cat COMMIT)" >&2
    mkdir -p ref
    ln -s "$(readlink -f ~{ref_fasta})" ref/ref.fa
    ln -s "$(readlink -f ~{ref_fasta_index})" ref/ref.fa.fai
    # Faithful port of the active (uncommented) body of nextflow/modules/introme_functions.nf.
    # MNV.sh scoring stays commented out upstream.
    echo "Current directory: $(pwd)"
    echo executing ls
    ls
    python3 AG_check/AG_check.py ~{variant_info_stripped} ref/ref.fa introme_annotate.vcf ~{prefix}.introme_annotate.ag_check.vcf
    bgzip -k ~{prefix}.introme_annotate.ag_check.vcf
    tabix -p vcf ~{prefix}.introme_annotate.ag_check.vcf.gz
    python3 ESE/scoring.py ~{variant_info_stripped} ./~{prefix}.introme_annotate.ESE.vcf ref/ref.fa
    mv ~{prefix}.introme_annotate.ESE.vcf ~{prefix}.introme_annotate.ESE.tsv
    bgzip -k ~{prefix}.introme_annotate.ESE.tsv
    # Upstream line is `tabix -s1 -b2 -e2 -c # file.tsv.gz`. In bash, # starts a
    # comment, so the filename never reaches tabix. -c '#' matches the #CHROM header.
    tabix -s1 -b2 -e2 -c '#' ~{prefix}.introme_annotate.ESE.tsv.gz
  >>>
  output {
    File ag_check = "~{prefix}.introme_annotate.ag_check.vcf.gz"
    File ag_check_index = "~{prefix}.introme_annotate.ag_check.vcf.gz.tbi"
    File ese_score = "~{prefix}.introme_annotate.ESE.tsv.gz"
    File ese_score_index = "~{prefix}.introme_annotate.ESE.tsv.gz.tbi"
  }
  runtime {
    docker: docker
    preemptible: preemptible
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: cpu
    bootDiskSizeGb: 30
  }
}

task SplicingAnno {
  input {
    File vcf
    File vcf_index
    File assets_tar
    File spliceai_vcf
    File mmsplice_vcf
    File pangolin_vcf
    File spip_vcf
    File spliceogen_tsv
    File ag_check
    File ag_check_index
    File ese_score
    File ese_score_index
    String prefix = "cohort"
    String genome_build = "hg38"
    String force_rerun_token = ""
    String docker = "headoncollusion/ag_check:v1.6"
    Int preemptible = 2
    Int cpu = 4
    Int memory_gb = 8
  }
  Int disk_gb = ceil(size(vcf, "GiB") + size(spliceai_vcf, "GiB") + size(mmsplice_vcf, "GiB") + size(pangolin_vcf, "GiB") + size(spip_vcf, "GiB") + size(spliceogen_tsv, "GiB") + size(assets_tar, "GiB") + 20)
  command <<<
    set -euo pipefail
    echo "force_rerun_token=~{force_rerun_token}" >&2
    tar -xzf ~{assets_tar}
    echo "introme_commit=$(cat COMMIT)" >&2
    # Faithful port of nextflow/modules/splicing_anno.nf.
    # Tool VCFs are localized inputs. Names below match vcfanno_splicing_ensemble.toml.
    # annotate.~{genome_build}.toml expects regions/, branchpointer/, and U12/ in -base-path.
    annotate_toml="tomls/annotate.~{genome_build}.toml"
    if [[ ! -f "$annotate_toml" ]]; then
      echo "No $annotate_toml in the pinned assets. This port stages hg38 only." >&2
      exit 1
    fi
    echo splicing_anno
    bgzip -c ~{spliceai_vcf} > spliceai.vcf.gz
    tabix -p vcf spliceai.vcf.gz
    bgzip -c ~{mmsplice_vcf} > mmsplice.vcf.gz
    tabix -p vcf mmsplice.vcf.gz
    bgzip -c ~{pangolin_vcf} > pangolin.vcf.gz
    tabix -p vcf pangolin.vcf.gz
    cp -L ~{spip_vcf} spip.in.vcf
    sed -i -E 's_(##SPiP output) (v[0-9]+.[0-9]*)_\1=\2_' spip.in.vcf
    bgzip -c spip.in.vcf > spip.vcf.gz
    tabix -p vcf spip.vcf.gz
    bgzip -c ~{spliceogen_tsv} > spliceogen.tsv.gz
    tabix -s1 -b2 -e3 spliceogen.tsv.gz
    ln -s ~{ag_check} introme_annotate.ag_check.vcf.gz
    ln -s ~{ag_check_index} introme_annotate.ag_check.vcf.gz.tbi
    ln -s ~{ese_score} introme_annotate.ESE.tsv.gz
    ln -s ~{ese_score_index} introme_annotate.ESE.tsv.gz.tbi
    pwd
    ls
    vcfanno \
        -base-path ./ \
        -p $(getconf _NPROCESSORS_ONLN) \
        -lua conf_ensemble.lua \
        "$annotate_toml" \
        ~{vcf} > step1.vcf
    bgzip step1.vcf
    vcfanno \
        -base-path ./ \
        -p $(getconf _NPROCESSORS_ONLN) \
        -lua conf_ensemble.lua \
        tomls/vcfanno_splicing_ensemble.toml \
        step1.vcf.gz > ~{prefix}.highquality.annotated.filtered.ensemblescored.vcf
    bgzip -k ~{prefix}.highquality.annotated.filtered.ensemblescored.vcf
  >>>
  output {
    File splicing_anno_vcf = "~{prefix}.highquality.annotated.filtered.ensemblescored.vcf.gz"
  }
  runtime {
    docker: docker
    preemptible: preemptible
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: cpu
    bootDiskSizeGb: 30
  }
}

task EnsembleInfer {
  input {
    File splicing_anno_vcf
    File assets_tar
    String prefix = "cohort"
    String commit_sha
    String force_rerun_token = ""
    String docker = "headoncollusion/ag_check:v1.6"
    Int preemptible = 1
    Int cpu = 4
    Int memory_gb = 16
  }
  Int disk_gb = ceil(size(splicing_anno_vcf, "GiB") + size(assets_tar, "GiB") + 20)
  command <<<
    set -euo pipefail
    echo "force_rerun_token=~{force_rerun_token}" >&2
    tar -xzf ~{assets_tar}
    echo "introme_commit=$(cat COMMIT)" >&2
    if [[ "$(tr -d '[:space:]' < COMMIT)" != "~{commit_sha}" ]]; then
      echo "assets COMMIT does not match fetch output ~{commit_sha}" >&2
      exit 1
    fi
    # Faithful port of nextflow/modules/ensemble.nf process ensemble_infer.
    # main.py imports train.py even in infer mode, so train.py is in the assets tar.
    # ensemble_train itself is not called.
    python3 ml/main.py infer \
        --model-path models/RandomForest_all_info_introme_train_run_model.pkl \
        --columns-path models/RandomForest_all_info_introme_train_run_columns.json \
        --input-vcf ~{splicing_anno_vcf} \
        --output-tsv ~{prefix}.introme.predictions.tsv
    bgzip -c ~{prefix}.introme.predictions.tsv > introme.predictions.tsv.gz
    tabix -S1 -s1 -b2 -e2 introme.predictions.tsv.gz
    introme_col=$(
        zcat introme.predictions.tsv.gz \
        | head -1 \
        | tr '\t' '\n' \
        | awk '$0 == "introme_score" {print NR; exit}'
    )
    sed "s/__INTROME_COL__/${introme_col}/" tomls/introme.toml > introme.generated.toml
    vcfanno \
        -base-path ./ \
        -p $(getconf _NPROCESSORS_ONLN) \
        introme.generated.toml \
        ~{splicing_anno_vcf} > ~{prefix}.introme.predictions.vcf
    bgzip ~{prefix}.introme.predictions.vcf
    tabix -p vcf ~{prefix}.introme.predictions.vcf.gz
    # pdvar second-hit lookup needs a gene and per-sample allele counts.
    # ensemble infer drops INFO/gene (commented out of INFO_FIELDS) and sets
    # format_fields=[], so the raw TSV has neither. Copy them from the annotated
    # VCF after vcfanno has already read the score column. Existing columns stay.
    export INTROME_TSV="~{prefix}.introme.predictions.tsv"
    export INTROME_VCF="~{prefix}.introme.predictions.vcf.gz"
    python3 -c 'import pathlib,sys,textwrap; pathlib.Path("enrich_carriers.py").write_text(textwrap.dedent(sys.stdin.read()))' <<'PY'
import csv
import gzip
import os
from pathlib import Path

tsv_path = Path(os.environ["INTROME_TSV"])
vcf_path = Path(os.environ["INTROME_VCF"])

def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open(encoding="utf-8", newline="")

def info_map(info):
    out = {}
    for part in info.split(";"):
        if not part or part == ".":
            continue
        if "=" in part:
            key, value = part.split("=", 1)
            out[key] = value
        else:
            out[part] = "1"
    return out

def allele_count(gt):
    gt = gt.split(":")[0].replace("|", "/")
    if gt in ("./.", ".|.", ".", ""):
        return 0
    alleles = [part for part in gt.split("/") if part not in ("", ".")]
    if not alleles or all(part == "0" for part in alleles):
        return 0
    return sum(1 for part in alleles if part != "0")

def norm_pos(value):
    text = str(value).strip()
    if text.endswith(".0"):
        text = text[:-2]
    return text

carriers = {}
genes = {}
samples = []
with open_text(vcf_path) as handle:
    for line in handle:
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            header = line.rstrip("\n").split("\t")
            samples = header[9:]
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 8:
            continue
        chrom, pos, ref, alt, info = cols[0], cols[1], cols[3], cols[4], cols[7]
        key = (chrom, norm_pos(pos), ref, alt.split(",")[0])
        gene = info_map(info).get("gene", "")
        symbols = [part for part in gene.replace(";", ",").split(",") if part and part != "."]
        genes[key] = symbols
        pairs = []
        if samples and len(cols) >= 10:
            fmt = cols[8].split(":")
            gt_index = fmt.index("GT") if "GT" in fmt else 0
            for sample, field in zip(samples, cols[9:]):
                bits = field.split(":")
                gt = bits[gt_index] if gt_index < len(bits) else field
                count = allele_count(gt)
                if count:
                    pairs.append((sample, count))
        carriers[key] = pairs

with tsv_path.open(encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    fieldnames = list(reader.fieldnames or [])
    rows = list(reader)
if "gene" not in fieldnames:
    fieldnames.append("gene")
if "samples_with_alt" not in fieldnames:
    fieldnames.append("samples_with_alt")
for row in rows:
    key = (row.get("CHROM", ""), norm_pos(row.get("POS", "")), row.get("REF", ""), row.get("ALT", ""))
    if not row.get("gene"):
        row["gene"] = repr(genes.get(key, []))
    if not row.get("samples_with_alt"):
        row["samples_with_alt"] = repr(carriers.get(key, []))
with tsv_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n", extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)
PY
    python3 enrich_carriers.py
    printf '%s\n' "~{commit_sha}" > introme.commit.txt
  >>>
  output {
    File introme_tsv = "~{prefix}.introme.predictions.tsv"
    File introme_vcf = "~{prefix}.introme.predictions.vcf.gz"
    File introme_vcf_index = "~{prefix}.introme.predictions.vcf.gz.tbi"
    File commit_txt = "introme.commit.txt"
  }
  runtime {
    docker: docker
    preemptible: preemptible
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: cpu
    bootDiskSizeGb: 30
  }
}
