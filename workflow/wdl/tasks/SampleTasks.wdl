version 1.0

# Tasks for one BAM (+ .bai). CRAM inputs are not supported.

task ExtractBamHeader {
  input {
    File bam
    File bam_index
    String sample_id
    String samtools_docker = "quay.io/biocontainers/samtools:1.24--h9dcdb79_1"
    Int preemptible = 3
  }
  Int disk_gb = ceil(size(bam, "GiB") * 2 + size(bam_index, "GiB") + 100)
  command <<<
    set -euo pipefail
    samtools view -H ~{bam} > ~{sample_id}.bam.header.txt
  >>>
  output {
    File header = "~{sample_id}.bam.header.txt"
  }
  runtime {
    docker: samtools_docker
    preemptible: preemptible
    memory: "4 GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: 1
  }
}

task CopyReferenceDict {
  input {
    File ref_dict
    String samtools_docker = "quay.io/biocontainers/samtools:1.24--h9dcdb79_1"
    Int preemptible = 3
  }
  command <<<
    set -euo pipefail
    cp "~{ref_dict}" reference.dict
  >>>
  output {
    File dict = "reference.dict"
  }
  runtime {
    docker: samtools_docker
    preemptible: preemptible
    memory: "4 GiB"
    disks: "local-disk 10 HDD"
    cpu: 1
  }
}

task SamtoolsReferenceDict {
  input {
    File ref_fasta
    String samtools_docker = "quay.io/biocontainers/samtools:1.24--h9dcdb79_1"
    Int preemptible = 3
  }
  Int disk_gb = ceil(size(ref_fasta, "GiB") * 2 + 100)
  command <<<
    set -euo pipefail
    samtools dict ~{ref_fasta} > reference.dict
  >>>
  output {
    File dict = "reference.dict"
  }
  runtime {
    docker: samtools_docker
    preemptible: preemptible
    memory: "4 GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: 1
  }
}

task FetchRepoAssets {
  input {
    String git_commit = "main"
    String python_docker = "python:3.11.11-slim-bookworm"
  }
  command <<<
    set -euo pipefail
    python3 - <<'PY'
import urllib.request

base = "https://raw.githubusercontent.com/gqworks/pdsample-terra/~{git_commit}"
pairs = [
    ("workflow/scripts/check_reference.py", "check_reference.py"),
    ("workflow/scripts/qc_summary.py", "qc_summary.py"),
    ("workflow/scripts/annotate_avi.py", "annotate_avi.py"),
    ("workflow/wdl/vendor/gatk-sv/format_svtk_vcf_for_gatk.py", "format_svtk_vcf_for_gatk.py"),
    ("workflow/assets/monopdaus-genes.grch38.bed", "panel.bed"),
    ("resources/panel/monopdaus-gene-inheritance.csv", "panel.csv"),
    ("resources/panel/monopdaus-genes.grch38.bed", "genes.bed"),
    ("resources/panel/monopdaus-exons.grch38.bed", "exons.bed"),
    ("workflow/assets/expansionhunter_pd.hg38.json", "eh_catalog.json"),
    ("config/config.yaml", "config.yaml"),
    ("config/columns.yaml", "columns.yaml"),
    ("config/gene_phenotypes.yaml", "phenotypes.yaml"),
    ("workflow/assets/verifybamid/1000g.phase3.100k.b38.vcf.gz.dat.UD", "svd.UD"),
    ("workflow/assets/verifybamid/1000g.phase3.100k.b38.vcf.gz.dat.V", "svd.V"),
    ("workflow/assets/verifybamid/1000g.phase3.100k.b38.vcf.gz.dat.bed", "svd.bed"),
    ("workflow/assets/verifybamid/1000g.phase3.100k.b38.vcf.gz.dat.mu", "svd.mu"),
]
for rel, dest in pairs:
    urllib.request.urlretrieve(base + "/" + rel, dest)
PY
  >>>
  output {
    File check_reference_py = "check_reference.py"
    File qc_summary_py = "qc_summary.py"
    File annotate_avi_py = "annotate_avi.py"
    File svtk_to_gatk_script = "format_svtk_vcf_for_gatk.py"
    File panel_bed = "panel.bed"
    File panel_csv = "panel.csv"
    File gene_bed = "genes.bed"
    File exon_bed = "exons.bed"
    File eh_catalog = "eh_catalog.json"
    File config_yaml = "config.yaml"
    File columns_yaml = "columns.yaml"
    File phenotypes_yaml = "phenotypes.yaml"
    File verifybamid_ud = "svd.UD"
    File verifybamid_v = "svd.V"
    File verifybamid_bed = "svd.bed"
    File verifybamid_mu = "svd.mu"
  }
  runtime {
    docker: python_docker
    memory: "2 GiB"
    disks: "local-disk 20 HDD"
    cpu: 1
  }
}

task CheckReference {
  input {
    File bam_header
    File reference_dict
    File check_reference_py
    String sample_id
    String python_docker = "python:3.11.11-slim-bookworm"
  }
  command <<<
    set -euo pipefail
    python3 ~{check_reference_py} --bam-header ~{bam_header} --reference-dict ~{reference_dict} --out ~{sample_id}.reference.json
  >>>
  output {
    File report = "~{sample_id}.reference.json"
  }
  runtime {
    docker: python_docker
    memory: "2 GiB"
    disks: "local-disk 10 HDD"
    cpu: 1
  }
}

task VerifyBamID {
  input {
    File bam
    File bam_index
    File ref_fasta
    File ref_fasta_index
    File svd_ud
    File svd_v
    File svd_bed
    File svd_mu
    String sample_id
    File reference_report
    String verifybamid_docker = "quay.io/biocontainers/verifybamid2:2.0.3--hc004090_0"
    Int preemptible = 3
  }
  # Cromwell copies the whole BAM and the reference onto this disk.
  Int disk_gb = ceil(size(bam, "GiB") * 2 + size(bam_index, "GiB") + size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + 200)
  command <<<
    set -euo pipefail
    align="~{bam}"
    index="~{bam_index}"
    if [[ "$align" != gs://* && -e "$index" ]]; then
      [[ -e "${align}.bai" ]] || ln -s "$index" "${align}.bai"
      [[ -e "${align%.bam}.bai" ]] || ln -s "$index" "${align%.bam}.bai"
    fi
    wc -c ~{reference_report} >/dev/null
    mkdir -p svd
    cp "~{svd_ud}" svd/resource.UD
    cp "~{svd_v}" svd/resource.V
    cp "~{svd_bed}" svd/resource.bed
    cp "~{svd_mu}" svd/resource.mu
    verifybamid2 --BamFile ~{bam} --Reference ~{ref_fasta} --SVDPrefix svd/resource --Output ~{sample_id}
  >>>
  output {
    File self_sm = "~{sample_id}.selfSM"
  }
  runtime {
    docker: verifybamid_docker
    preemptible: preemptible
    memory: "8 GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: 2
  }
}

task SomalierCheck {
  input {
    File bam
    File bam_index
    File ref_fasta
    File ref_fasta_index
    File sites
    String sample_id
    File reference_report
    String somalier_docker = "quay.io/biocontainers/somalier:0.3.5--h5205c93_0"
    Int preemptible = 3
  }
  Int disk_gb = ceil(size(bam, "GiB") * 2 + size(bam_index, "GiB") + size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + 200)
  command <<<
    set -euo pipefail
    align="~{bam}"
    index="~{bam_index}"
    if [[ "$align" != gs://* && -e "$index" ]]; then
      [[ -e "${align}.bai" ]] || ln -s "$index" "${align}.bai"
      [[ -e "${align%.bam}.bai" ]] || ln -s "$index" "${align%.bam}.bai"
    fi
    wc -c ~{reference_report} >/dev/null
    mkdir -p extract
    somalier extract -d extract --sites ~{sites} -f ~{ref_fasta} ~{bam}
    extract_file=$(find extract -name '*.somalier' | head -1)
    test -n "$extract_file"
    # One extract has no pairs. somalier relate still writes the sample's sex fields.
    somalier relate "$extract_file"
    if [[ ! -f relatedness.samples.tsv ]]; then
      found=$(find . -name '*samples*.tsv' | head -1)
      test -n "$found"
      cp "$found" relatedness.samples.tsv
    fi
    cp relatedness.samples.tsv ~{sample_id}.somalier.samples.tsv
  >>>
  output {
    File samples_tsv = "~{sample_id}.somalier.samples.tsv"
  }
  runtime {
    docker: somalier_docker
    preemptible: preemptible
    memory: "4 GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: 1
  }
}

task Mosdepth {
  input {
    File bam
    File bam_index
    File ref_fasta
    File ref_fasta_index
    File bed
    String sample_id
    File reference_report
    String mosdepth_docker = "quay.io/biocontainers/mosdepth:0.3.14--h87be163_2"
    Int preemptible = 3
    Int threads = 4
  }
  Int disk_gb = ceil(size(bam, "GiB") * 2 + size(bam_index, "GiB") + size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + 200)
  command <<<
    set -euo pipefail
    align="~{bam}"
    index="~{bam_index}"
    if [[ "$align" != gs://* && -e "$index" ]]; then
      [[ -e "${align}.bai" ]] || ln -s "$index" "${align}.bai"
      [[ -e "${align%.bam}.bai" ]] || ln -s "$index" "${align%.bam}.bai"
    fi
    wc -c ~{reference_report} >/dev/null
    mosdepth -t ~{threads} -n -f ~{ref_fasta} --thresholds 1,10,20,30 --by ~{bed} ~{sample_id} ~{bam}
  >>>
  output {
    File summary = "~{sample_id}.mosdepth.summary.txt"
  }
  runtime {
    docker: mosdepth_docker
    preemptible: preemptible
    memory: "8 GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: threads
  }
}

task QcGate {
  input {
    File qc_summary_py
    File selfsm
    File somalier_samples
    File mosdepth_summary
    File reference_report
    String sample_id
    String provided_sex = ""
    Float max_freemix = 0.02
    Float min_mean_coverage = 0
    Boolean hard_fail_contamination = true
    Boolean hard_fail_sex = false
    Boolean hard_fail_coverage = false
    String python_docker = "python:3.11.11-slim-bookworm"
  }
  command <<<
    set -euo pipefail
    python3 ~{qc_summary_py} \
      --sample-id ~{sample_id} \
      --selfsm ~{selfsm} \
      --somalier ~{somalier_samples} \
      --mosdepth ~{mosdepth_summary} \
      --reference ~{reference_report} \
      --provided-sex "~{provided_sex}" \
      --max-freemix ~{max_freemix} \
      --min-mean-coverage ~{min_mean_coverage} \
      --hard-fail-contamination ~{true="true" false="false" hard_fail_contamination} \
      --hard-fail-sex ~{true="true" false="false" hard_fail_sex} \
      --hard-fail-coverage ~{true="true" false="false" hard_fail_coverage} \
      --out qc.json \
      --sex-out sex.txt
    touch qc.ok
  >>>
  output {
    File qc_json = "qc.json"
    File gate = "qc.ok"
    String sex = read_lines("sex.txt")[0]
  }
  runtime {
    docker: python_docker
    memory: "2 GiB"
    disks: "local-disk 10 HDD"
    cpu: 1
  }
}

task HaplotypeCaller {
  input {
    File bam
    File bam_index
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File? interval_bed
    Boolean use_intervals = true
    String sample_id
    File qc_gate
    String? requester_pays_project
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
    Int preemptible = 3
    Int memory_mb = 32768
  }
  Int disk_gb = ceil(size(bam, "GiB") * 3 + size(bam_index, "GiB") + size(ref_fasta, "GiB") + size(ref_dict, "GiB") + 250)
  String rp = select_first([requester_pays_project, ""])
  command <<<
    set -euo pipefail
    interval_args=()
    if [[ "~{use_intervals}" == "true" ]]; then
      interval_args=(-L ~{interval_bed})
    fi
    wc -c ~{qc_gate} >/dev/null
    mem_unit=$(echo "${MEM_UNIT:-}" | cut -c 1)
    if [[ "$mem_unit" == "M" ]]; then
      available_memory_mb=$(awk "BEGIN {print int($MEM_SIZE)}")
    elif [[ "$mem_unit" == "G" ]]; then
      available_memory_mb=$(awk "BEGIN {print int($MEM_SIZE * 1024)}")
    else
      available_memory_mb=$(free -m | awk '/^Mem/ {print $2}')
    fi
    java_memory_size_mb=$((available_memory_mb - 1024))
    rp_arg=""
    if [[ -n "~{rp}" ]]; then
      rp_arg="--gcs-project-for-requester-pays ~{rp}"
    fi
    gatk --java-options "-Xmx${java_memory_size_mb}m -Xms${java_memory_size_mb}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
      HaplotypeCaller \
      $rp_arg \
      -R ~{ref_fasta} \
      -I ~{bam} \
      "${interval_args[@]}" \
      -O ~{sample_id}.g.vcf.gz \
      -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation \
      -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
      -ERC GVCF
  >>>
  output {
    File gvcf = "~{sample_id}.g.vcf.gz"
    File gvcf_index = "~{sample_id}.g.vcf.gz.tbi"
  }
  runtime {
    docker: gatk_docker
    preemptible: preemptible
    memory: "~{memory_mb} MiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: 4
    bootDiskSizeGb: 15
  }
  parameter_meta {
    bam: { localization_optional: true }
    bam_index: { localization_optional: true }
    ref_fasta: { localization_optional: true }
    ref_fasta_index: { localization_optional: true }
    ref_dict: { localization_optional: true }
  }
}

task GenotypeHardFilter {
  input {
    File gvcf
    File gvcf_index
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File? interval_bed
    Boolean use_intervals = true
    String sample_id
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
    Int preemptible = 3
  }
  Int disk_gb = ceil(size(gvcf, "GiB") * 4 + size(gvcf_index, "GiB") + size(ref_fasta, "GiB") + 200)
  command <<<
    set -euo pipefail
    interval_args=()
    if [[ "~{use_intervals}" == "true" ]]; then
      interval_args=(-L ~{interval_bed})
    fi
    # GATK requires the sequence dictionary beside the FASTA; Cromwell localizes dict separately.
    ref_bundle_dir="gatk_ref"
    mkdir -p "$ref_bundle_dir"
    ref_base=$(basename ~{ref_fasta})
    ln -sf ~{ref_fasta} "$ref_bundle_dir/$ref_base"
    ln -sf ~{ref_fasta_index} "$ref_bundle_dir/${ref_base}.fai"
    dict_name="${ref_base%.*}.dict"
    ln -sf ~{ref_dict} "$ref_bundle_dir/$dict_name"
    ref_fasta_path="$ref_bundle_dir/$ref_base"
    # GATK hard filters for a single sample. VQSR needs a cohort to fit
    # tranches, and this workflow does not ship a CNN score model.
    gatk --java-options "-Xmx6g" GenotypeGVCFs -R "$ref_fasta_path" -V ~{gvcf} "${interval_args[@]}" -O raw.vcf.gz
    gatk --java-options "-Xmx6g" SelectVariants -V raw.vcf.gz -select-type SNP -O snps.vcf.gz
    gatk --java-options "-Xmx6g" VariantFiltration -V snps.vcf.gz \
      --filter-expression "QD < 2.0 || QUAL < 30.0 || SOR > 3.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
      --filter-name "hard_filter" -O snps.filtered.vcf.gz
    gatk --java-options "-Xmx6g" SelectVariants -V raw.vcf.gz -select-type INDEL -O indels.vcf.gz
    gatk --java-options "-Xmx6g" VariantFiltration -V indels.vcf.gz \
      --filter-expression "QD < 2.0 || QUAL < 30.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
      --filter-name "hard_filter" -O indels.filtered.vcf.gz
    gatk --java-options "-Xmx4g" MergeVcfs -I snps.filtered.vcf.gz -I indels.filtered.vcf.gz -O ~{sample_id}.genotype.vcf.gz
  >>>
  output {
    File vcf = "~{sample_id}.genotype.vcf.gz"
    File vcf_index = "~{sample_id}.genotype.vcf.gz.tbi"
  }
  runtime {
    docker: gatk_docker
    preemptible: preemptible
    memory: "8 GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: 2
  }
}

task NormPass {
  input {
    File vcf
    File vcf_index
    File ref_fasta
    File ref_fasta_index
    String sample_id
    String bcftools_docker = "quay.io/biocontainers/bcftools:1.24--h118bc1c_4"
    Int preemptible = 3
  }
  Int disk_gb = ceil(size(vcf, "GiB") * 3 + size(vcf_index, "GiB") + size(ref_fasta, "GiB") + 150)
  command <<<
    set -euo pipefail
    bcftools view -f PASS ~{vcf} \
      | bcftools norm -m -both -f ~{ref_fasta} -Oz -o ~{sample_id}.small.vcf.gz
    bcftools index -t ~{sample_id}.small.vcf.gz
  >>>
  output {
    File normalized_vcf = "~{sample_id}.small.vcf.gz"
    File normalized_vcf_index = "~{sample_id}.small.vcf.gz.tbi"
  }
  runtime {
    docker: bcftools_docker
    preemptible: preemptible
    memory: "4 GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: 1
  }
}

task VepAnnotate {
  input {
    File vcf
    File vcf_index
    File ref_fasta
    File ref_fasta_index
    File cache_dir
    File? gnomad_vcf
    File? gnomad_vcf_index
    File? clinvar_vcf
    File? clinvar_vcf_index
    String extra_args = ""
    String vep_docker = "quay.io/biocontainers/ensembl-vep:116.2--pl5321h2a3209d_0"
    Int preemptible = 3
    Int threads = 4
  }
  Int disk_gb = ceil(size(vcf, "GiB") + size(vcf_index, "GiB") + size(ref_fasta, "GiB") + size(cache_dir, "GiB") * 4 + 200)
  command <<<
    set -euo pipefail
    custom_args=()
    if [[ -n "~{gnomad_vcf}" ]]; then
      test -s ~{gnomad_vcf_index}
      custom_args+=(--custom ~{gnomad_vcf},gnomAD,vcf,exact,0,AF,AF_popmax,nhomalt)
    fi
    if [[ -n "~{clinvar_vcf}" ]]; then
      test -s ~{clinvar_vcf_index}
      custom_args+=(--custom ~{clinvar_vcf},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT)
    fi
    cache_arg="~{cache_dir}"
    case "~{cache_dir}" in
      *.tar.gz|*.tgz)
        mkdir -p vep_cache
        tar -xzf "~{cache_dir}" -C vep_cache
        cache_arg="vep_cache"
        ;;
    esac
    vep --offline --cache --dir_cache "$cache_arg" --fasta ~{ref_fasta} --assembly GRCh38 \
      --vcf --allele_number --symbol --hgvs --canonical --pick --force_overwrite \
      --fork ~{threads} \
      "${custom_args[@]}" ~{extra_args} \
      -i ~{vcf} -o panel.vep.vcf
    bgzip -f panel.vep.vcf
    tabix -p vcf panel.vep.vcf.gz
  >>>
  output {
    File vep_vcf = "panel.vep.vcf.gz"
    File vep_index = "panel.vep.vcf.gz.tbi"
  }
  runtime {
    docker: vep_docker
    preemptible: preemptible
    memory: "16 GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: threads
  }
}

task LiftGnomadInfo {
  input {
    File vcf
    File vcf_index
    String bcftools_docker = "quay.io/biocontainers/bcftools:1.24--h118bc1c_4"
    Int preemptible = 3
  }
  Int disk_gb = ceil(size(vcf, "GiB") * 3 + size(vcf_index, "GiB") + 150)
  command <<<
    set -euo pipefail
    header=$(bcftools view -h ~{vcf})
    if grep -q '^##INFO=<ID=gnomAD_AF,' <<<"$header"; then
      bcftools view -Oz -o lifted.vcf.gz ~{vcf}
      bcftools index -t lifted.vcf.gz
      exit 0
    fi
    if ! grep -q '^##INFO=<ID=gnomAD,' <<<"$header"; then
      echo "No VEP custom gnomAD INFO field. Introme rare filter stays a no-op." >&2
      bcftools view -Oz -o lifted.vcf.gz ~{vcf}
      bcftools index -t lifted.vcf.gz
      exit 0
    fi
    bcftools view ~{vcf} | awk 'BEGIN { OFS="\t" }
      /^##/ { print; next }
      /^#CHROM/ {
        print "##INFO=<ID=gnomAD_AF,Number=1,Type=Float,Description=\"AF copied from VEP custom INFO/gnomAD field 1\">"
        print "##INFO=<ID=AF_popmax,Number=1,Type=Float,Description=\"AF_popmax copied from VEP custom INFO/gnomAD field 2\">"
        print
        next
      }
      {
        info = $8
        val = ""
        n = split(info, parts, ";")
        for (i = 1; i <= n; i++) {
          if (parts[i] ~ /^gnomAD=/) val = substr(parts[i], 8)
        }
        if (val != "" && val != ".") {
          split(val, fields, "|")
          af = fields[1]
          pop = fields[2]
          sub(/^[A-Za-z0-9_]+=/, "", af)
          sub(/^[A-Za-z0-9_]+=/, "", pop)
          if (af != "" && af != ".") info = info ";gnomAD_AF=" af
          if (pop != "" && pop != ".") info = info ";AF_popmax=" pop
          $8 = info
        }
        print
      }' | bcftools view -Oz -o lifted.vcf.gz
    bcftools index -t lifted.vcf.gz
  >>>
  output {
    File lifted_vcf = "lifted.vcf.gz"
    File lifted_vcf_index = "lifted.vcf.gz.tbi"
  }
  runtime {
    docker: bcftools_docker
    preemptible: preemptible
    memory: "4 GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: 1
  }
}

task ExpansionHunter {
  input {
    File bam
    File bam_index
    File ref_fasta
    File ref_fasta_index
    File catalog
    String sample_id
    String sex = ""
    File qc_gate
    String eh_docker = "quay.io/biocontainers/expansionhunter:5.0.0--hc26b3af_5"
    Int preemptible = 3
    Int threads = 4
  }
  Int disk_gb = ceil(size(bam, "GiB") * 2 + size(bam_index, "GiB") + size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + 200)
  command <<<
    set -euo pipefail
    align="~{bam}"
    index="~{bam_index}"
    if [[ "$align" != gs://* && -e "$index" ]]; then
      [[ -e "${align}.bai" ]] || ln -s "$index" "${align}.bai"
      [[ -e "${align%.bam}.bai" ]] || ln -s "$index" "${align%.bam}.bai"
    fi
    wc -c ~{qc_gate} >/dev/null
    sex_args=()
    if [[ "~{sex}" == "male" || "~{sex}" == "female" ]]; then
      sex_args=(--sex "~{sex}")
    fi
    ExpansionHunter \
      --reads ~{bam} \
      --reference ~{ref_fasta} \
      --variant-catalog ~{catalog} \
      --output-prefix ~{sample_id} \
      --threads ~{threads} \
      "${sex_args[@]}"
    bgzip -f ~{sample_id}.vcf
    tabix -p vcf ~{sample_id}.vcf.gz
  >>>
  output {
    File vcf = "~{sample_id}.vcf.gz"
    File vcf_index = "~{sample_id}.vcf.gz.tbi"
    File json = "~{sample_id}.json"
  }
  runtime {
    docker: eh_docker
    preemptible: preemptible
    memory: "16 GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: threads
  }
}

task Gauchian {
  input {
    File bam
    File bam_index
    File ref_fasta
    File ref_fasta_index
    String sample_id
    File qc_gate
    String gauchian_docker = "python:3.11.11-slim-bookworm"
    Int preemptible = 3
    Int threads = 4
  }
  Int disk_gb = ceil(size(bam, "GiB") * 2 + size(bam_index, "GiB") + size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + 200)
  command <<<
    set -euo pipefail
    align="~{bam}"
    index="~{bam_index}"
    if [[ "$align" != gs://* && -e "$index" ]]; then
      [[ -e "${align}.bai" ]] || ln -s "$index" "${align}.bai"
      [[ -e "${align%.bam}.bai" ]] || ln -s "$index" "${align%.bam}.bai"
    fi
    wc -c ~{qc_gate} >/dev/null
    if ! command -v gauchian >/dev/null; then
      pip install --no-cache-dir gauchian==1.0.2
    fi
    echo "~{bam}" > manifest.txt
    mkdir -p out
    gauchian --manifest manifest.txt --genome 38 --prefix ~{sample_id} --outDir out --threads ~{threads} --reference ~{ref_fasta}
    json=$(find out -name '*.json' -print -quit)
    if [[ -z "$json" ]]; then
      echo "Gauchian wrote no JSON under out/" >&2
      find out -type f >&2 || true
      exit 1
    fi
    cp "$json" ~{sample_id}.gauchian.json
  >>>
  output {
    File json = "~{sample_id}.gauchian.json"
  }
  runtime {
    docker: gauchian_docker
    preemptible: preemptible
    memory: "8 GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: threads
  }
}

task ReleaseAfterQc {
  input {
    File payload
    File gate
    String name
  }
  Int disk_gb = ceil(size(payload, "GiB") * 2 + 200)
  command <<<
    set -euo pipefail
    wc -c ~{gate} >/dev/null
    ln ~{payload} ~{name} || cp -a ~{payload} ~{name}
  >>>
  output {
    File released = "~{name}"
  }
  runtime {
    docker: "python:3.11.11-slim-bookworm"
    memory: "2 GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: 1
  }
}

task AnnotateAvi {
  input {
    File vcf
    File annotate_avi_py
    # AlphaGenome "AVI SNV scores" download (bgzipped TSV) and its Tabix index.
    # Unset writes a header-only table and leaves AVI not_assessed.
    File? avi_table
    File? avi_table_index
    String python_docker = "python:3.11.11-slim-bookworm"
    Int preemptible = 3
  }
  Int disk_gb = ceil(size(vcf, "GiB") + size(avi_table, "GiB") + size(avi_table_index, "GiB") + 150)
  command <<<
    set -euo pipefail
    # Local Tabix lookup. No API key and no AlphaGenome network call.
    table='~{default="" avi_table}'
    index='~{default="" avi_table_index}'
    args=()
    if [[ -n "${table}" ]]; then
      pip install --no-cache-dir "pysam==0.24.1"
      args+=(--table "${table}")
      # Cromwell may localize the index in a different directory than the table.
      if [[ -n "${index}" ]]; then
        args+=(--index "${index}")
      fi
    fi
    python3 "~{annotate_avi_py}" --vcf "~{vcf}" --out avi.tsv ${args[@]+"${args[@]}"}
  >>>
  output {
    File avi_tsv = "avi.tsv"
  }
  runtime {
    docker: python_docker
    preemptible: preemptible
    memory: "8 GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: 2
  }
}

task Assess {
  input {
    String sample_id
    String sex
    File vep_vcf
    File small_vcf
    File sv_vcf
    File expansionhunter_vcf
    File gauchian_json
    File mtdna_vcf
    File introme_tsv
    File? avi_tsv
    File qc_json
    File config_yaml
    File columns_yaml
    File phenotypes_yaml
    File panel_csv
    File? gene_bed
    File? exon_bed
    File? expression_tsv
    File? splicevault_tsv
    File? prs_weights
    String? age_onset
    String? symptoms
    String? ethnicity
    File? phenotype_report
    String pdvar_docker = "python:3.11.11-slim-bookworm"
    String pdvar_git_url = "https://github.com/gqworks/pdsample-terra.git"
    String pdvar_git_commit = "main"
  }
  String age_value = select_first([age_onset, ""])
  String symptoms_value = select_first([symptoms, ""])
  String ethnicity_value = select_first([ethnicity, ""])
  # write_json keeps quotes and newlines in free-text symptoms out of the shell.
  File phenotype_json = write_json({
    "age_onset": age_value,
    "symptoms": symptoms_value,
    "ethnicity": ethnicity_value
  })
  Int disk_gb = ceil(size(vep_vcf, "GiB") + size(small_vcf, "GiB") + size(sv_vcf, "GiB") + size(expansionhunter_vcf, "GiB") + size(mtdna_vcf, "GiB") + size(introme_tsv, "GiB") + size(select_first([avi_tsv, vep_vcf]), "GiB") + 150)
  command <<<
    set -euo pipefail
    if ! python3 -c "import pdvar" >/dev/null 2>&1; then
      if [[ -z "~{pdvar_git_commit}" ]]; then
        echo "pdvar is not in the image. Set pdvar_git_commit or pdvar_docker to an image built from workflow/containers/pdvar." >&2
        exit 1
      fi
      apt-get update
      apt-get install -y --no-install-recommends git ca-certificates
      rm -rf /var/lib/apt/lists/*
      pip install --no-cache-dir "git+~{pdvar_git_url}@~{pdvar_git_commit}"
    fi
    export SAMPLE_ID="~{sample_id}"
    export SEX="~{sex}"
    export VEP_VCF="~{vep_vcf}"
    export SMALL_VCF="~{small_vcf}"
    export SV_VCF="~{sv_vcf}"
    export EH_VCF="~{expansionhunter_vcf}"
    export GAUCHIAN_JSON="~{gauchian_json}"
    export MTDNA_VCF="~{mtdna_vcf}"
    export INTROME_TSV="~{introme_tsv}"
    export QC_JSON="~{qc_json}"
    export GENE_BED="~{gene_bed}"
    export EXON_BED="~{exon_bed}"
    export EXPRESSION="~{expression_tsv}"
    export SPLICEVAULT="~{splicevault_tsv}"
    export PRS_WEIGHTS="~{prs_weights}"
    export AVI_TSV="~{avi_tsv}"
    export PHENOTYPE_JSON="~{phenotype_json}"
    export PHENOTYPE_REPORT="~{phenotype_report}"
    python3 - <<'PY'
import json, os
from pathlib import Path

def maybe(key):
    value = os.environ.get(key, "").strip()
    return value or None

qc = json.loads(Path(os.environ["QC_JSON"]).read_text(encoding="utf-8"))
raw_phenotype = json.loads(Path(os.environ["PHENOTYPE_JSON"]).read_text(encoding="utf-8"))
phenotype = {}
for key in ("age_onset", "symptoms", "ethnicity"):
    value = str(raw_phenotype.get(key) or "").strip()
    if value:
        phenotype[key] = value
report_path = os.environ.get("PHENOTYPE_REPORT", "").strip()
if report_path:
    phenotype["report_text"] = Path(report_path).read_text(encoding="utf-8")
document = {
    "sample_id": os.environ["SAMPLE_ID"],
    "sex": os.environ.get("SEX") or qc.get("sex", {}).get("value") or "unknown",
    "vep_vcf": os.environ["VEP_VCF"],
    "small_vcf": os.environ["SMALL_VCF"],
    "sv_vcf": os.environ["SV_VCF"],
    "sv_status": "called",
    "expansionhunter_vcf": os.environ["EH_VCF"],
    "gauchian_json": os.environ["GAUCHIAN_JSON"],
    "mtdna_vcf": os.environ["MTDNA_VCF"],
    "introme_tsv": os.environ["INTROME_TSV"],
    "gene_bed": maybe("GENE_BED"),
    "exon_bed": maybe("EXON_BED"),
    "expression": maybe("EXPRESSION"),
    "splicevault": maybe("SPLICEVAULT"),
    "prs_weights": maybe("PRS_WEIGHTS"),
    "avi_tsv": maybe("AVI_TSV"),
    "phenotype": phenotype,
    "qc": qc,
}
Path("sample.json").write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
PY
    extra=()
    if [[ -n "$GENE_BED" ]]; then extra+=(--gene-bed "$GENE_BED"); fi
    if [[ -n "$EXON_BED" ]]; then extra+=(--exon-bed "$EXON_BED"); fi
    pdvar assess sample.json \
      --config ~{config_yaml} \
      --columns ~{columns_yaml} \
      --phenotypes ~{phenotypes_yaml} \
      --panel ~{panel_csv} \
      --outdir out \
      "${extra[@]}"
  >>>
  output {
    File sample_json = "sample.json"
    File candidates_tsv = "out/candidates.tsv"
    File report_html = "out/report.html"
    File report_md = "out/report.md"
    File qc_summary = "out/qc_summary.json"
    File liability_tsv = "out/liability.tsv"
    File followup_tsv = "out/followup.tsv"
  }
  runtime {
    docker: pdvar_docker
    memory: "8 GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: 2
  }
}
