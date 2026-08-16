#!/usr/bin/env python3
"""Build reader-facing supplementary methods tables from audited TSV sources.

Outputs:
  Supplementary Table 1: original dataset summary
  Supplementary Table 2: original reproductive-age-bin definitions and counts
  Supplementary Table 3a: 117 frozen model-development samples
  Supplementary Table 3b: 300 external-projection samples
  Supplementary Table 4: external-projection contrast definitions
  Supplementary Table 5: integrated-clock results for 34 primary contrasts

Only Python's standard library is required. The XLSX writer deliberately uses
inline strings so the workbook is reproducible without optional spreadsheet
packages.
"""

from __future__ import annotations

import csv
import math
import re
import zipfile
from datetime import date
from pathlib import Path
from xml.sax.saxutils import escape


PROJECT = Path(__file__).resolve().parents[1]
METHODS = PROJECT / "methods_tables"
OUT = METHODS / "publication_supplementary_tables"


REFERENCES = {
    "CRA003645": (
        "Zhou et al., 2021",
        "Temporal transcriptomic landscape of postnatal mouse ovaries reveals dynamic gene signatures associated with ovarian aging",
        "10.1093/hmg/ddab163",
        "publication_verified",
    ),
    "CRA017492": (
        "Liang et al., 2024",
        "Nicotinamide mononucleotide supplementation rescues mitochondrial and energy metabolism functions and ameliorates inflammatory states in the ovaries of aging mice",
        "10.1002/mco2.727",
        "publication_verified",
    ),
    "GSE106339": (
        "Sominsky et al., 2019",
        "Acylated Ghrelin Supports the Ovarian Transcriptome and Follicles in the Mouse: Implications for Fertility",
        "10.3389/fendo.2018.00815",
        "publication_verified",
    ),
    "GSE128240": (
        "Li et al., 2019",
        "Transcriptome analysis reveals that cyclophosphamide induces premature ovarian failure by blocking cholesterol biosynthesis pathway",
        "10.1016/j.lfs.2019.116999",
        "publication_verified",
    ),
    "GSE154890": (
        "Ansere et al., 2021",
        "Cellular hallmarks of aging emerge in the ovary prior to primordial follicle depletion",
        "10.1016/j.mad.2020.111425",
        "publication_verified",
    ),
    "GSE194196": (
        "Wang et al., 2022",
        "BNC1 deficiency-triggered ferroptosis through the NF2-YAP pathway induces primary ovarian insufficiency",
        "10.1038/s41467-022-33323-8",
        "publication_verified",
    ),
    "GSE228975": (
        "Yao et al., 2023",
        "SRSF1 is essential for primary follicle development by regulating granulosa cell survival via mRNA alternative splicing",
        "10.1007/s00018-023-04979-2",
        "publication_verified",
    ),
    "GSE235021": (
        "Tu et al., 2026",
        "Whole-grain black rice diet rewires the single-cell transcriptomic landscape of age-related ovarian decline in mice",
        "10.1007/s11427-025-3046-4",
        "publication_verified",
    ),
    "GSE240547": (
        "Liu et al., 2024",
        "High-fat and high-sucrose diet impairs female reproduction by altering ovarian transcriptomic and metabolic signatures",
        "10.1186/s12967-024-04952-y",
        "publication_verified",
    ),
    "GSE249151": (
        "Li et al., 2024",
        "Transcriptome analysis during 4-vinylcyclohexene diepoxide exposure-induced premature ovarian insufficiency in mice",
        "10.7717/peerj.17251",
        "publication_verified",
    ),
    "GSE252791": (
        "Song et al., 2024",
        "Deleterious variants in RNF111 impair female fertility and induce premature ovarian insufficiency in humans and mice",
        "10.1007/s11427-024-2606-6",
        "publication_verified",
    ),
    "GSE273724": (
        "Gui et al., 2025",
        "HEP14 treatment improves ovarian function in aged mice through mitophagy enhancement and oxidative stress reduction",
        "10.1038/s42003-025-08576-w",
        "publication_verified",
    ),
    "GSE279446": (
        "Stefano et al., 2026",
        "Derepression of a single microRNA target causes female infertility in mice",
        "10.1093/nar/gkaf1357",
        "publication_verified",
    ),
    "GSE279573": (
        "Yamanashi et al., 2024",
        "Dietary phytosterols induce infertility in female mice via epigenomic modulations",
        "10.1038/s42003-024-07233-y",
        "publication_verified",
    ),
    "GSE329024": (
        "Yang et al., 2026",
        "Distinct HPO axis responses and ovarian aging trajectories to chronic unpredictable mild stress in reproductively young versus middle-aged female mice",
        "10.64898/2026.04.24.720585",
        "preprint_verified",
    ),
    "GSE98497": (
        "Guo et al., 2018",
        "Oocyte stage-specific effects of MTOR determine granulosa cell fate and oocyte quality in mice",
        "10.1073/pnas.1800352115",
        "publication_verified",
    ),
    "NIHMS1554107": (
        "Zhang et al., 2020",
        "Inflammaging is associated with shifted macrophage ontogeny and polarization in the aging mouse ovary",
        "10.1530/REP-19-0330",
        "publication_verified",
    ),
    "PRJNA1076509": (
        "Kim et al., 2026",
        "Estropausal gut microbiota transplant improves measures of ovarian function in adult mice",
        "10.1038/s43587-026-01069-3",
        "publication_verified",
    ),
    "PRJNA787028": (
        "Kim and You, 2022",
        "High Housing Density-Induced Chronic Stress Diminishes Ovarian Reserve via Granulosa Cell Apoptosis by Angiotensin II Overexpression in Mice",
        "10.3390/ijms23158614",
        "publication_verified",
    ),
    "PRJNA924255": (
        "Chen et al., 2023",
        "Semaphorin 4C regulates ovarian steroidogenesis through RHOA/ROCK1-mediated actin cytoskeleton rearrangement",
        "10.1093/molehr/gaad010",
        "publication_verified",
    ),
    "PRJNA1002222": (
        "Yang et al., 2024",
        "NADase CD38 is a key determinant of ovarian aging",
        "10.1038/s43587-023-00532-9",
        "publication_verified",
    ),
    "2022_eLife_Morris_cycling_murine_ovary_atlas": (
        "Morris et al., 2022",
        "A single-cell atlas of the cycling murine ovary",
        "10.7554/eLife.77239",
        "publication_verified",
    ),
    "2024_Cell_Winkler_FRT_scRNA_spatial": (
        "Winkler et al., 2024",
        "The cycling and aging mouse female reproductive tract at single-cell resolution",
        "10.1016/j.cell.2024.01.021",
        "publication_verified",
    ),
    "2024_LifeMedicine_Jiang_whole_life_mouse_ovary": (
        "Jiang et al., 2024",
        "Single-cell RNA sequencing reveals the intercellular crosstalk and the regulatory landscape of stromal cells during the whole life of the mouse ovary",
        "10.1093/lifemedi/lnae041",
        "publication_verified",
    ),
    "2024_NatureAging_Isola_aging_mouse_ovary_atlas": (
        "Isola et al., 2024",
        "A single-cell atlas of the aging mouse ovary",
        "10.1038/s43587-023-00552-5",
        "publication_verified",
    ),
    "2025_CommBio_Zhao_developing_ovary_scRNA_spatial": (
        "Zhao et al., 2025",
        "Spatiotemporal and single-cell atlases to dissect regional specific cell types of the developing ovary",
        "10.1038/s42003-025-08277-4",
        "publication_verified",
    ),
    "2025_GeroScience_Zhang_multiomic_aging_mouse_ovary": (
        "Zhang et al., 2025",
        "A multi-omic single-cell landscape of the aging mouse ovary",
        "10.1007/s11357-025-01556-2",
        "publication_verified",
    ),
    "2026_CommBio_Galligos_ovarian_inflammaging": (
        "Galligos et al., 2026",
        "Multicellular origins of murine ovarian inflammaging",
        "10.1038/s42003-026-09826-1",
        "publication_verified",
    ),
    "PRJNA863443": (
        "Kim et al., 2025",
        "Systematic characterization of the ovarian landscape across mouse menopause models",
        "10.1101/2025.06.10.658920",
        "preprint_verified",
    ),
}


LIMITATIONS = {
    "CRA003645": "Neonatal and prepubertal samples are developmental/OOD; adult samples provide the held-out age-placement comparison.",
    "GSE106339": "GOAT-knockout samples are classified as genetic models and WT samples as reference controls within the analysis-domain sample set; all are prepubertal (3 weeks), so the contrast remains descriptive and excluded from primary adult contrast analyses.",
    "GSE154890": "Natural-age placement cohort without a treatment contrast.",
    "GSE228975": "Neonatal granulosa-cell knockout cohort; retained as developmental/OOD descriptive projection.",
    "GSE235021": "Diet effects are matched within age using polished-rice controls pooled across bulk RNA-seq and scRNA-derived pseudobulk; 6- and 9-month groups are not pooled, and mixed-modality effects are modality/batch-confounded by design.",
    "GSE249151": "The day-30 VCD comparison lacks an age-matched day-30 vehicle group and is flagged as control-age-confounded.",
    "GSE273724": "Natural-age placement and treatment comparisons are reported separately.",
    "GSE279573": "The KO oil-free versus KO chow contrast is diet/fertility-confounded.",
    "GSE329024": "Preprint cohort; age-stratified stress comparisons are retained separately.",
    "PRJNA863443": "scRNA-derived pseudobulk used for external projection only; several matched-control groups contain two or three animals.",
    "PRJNA924255": "Curated in the metadata inventory but not included in the current external-projection contrast manifest.",
    "PRJNA1002222": "Cd38-knockout ovaries are compared with same-age WT controls separately at 8 and 12 months; the 2-month WT group is a young reference and is not included in either knockout contrast.",
}


PRIMARY_CLASSES = {
    "intervention",
    "diet",
    "genetic_model",
    "stress_model",
    "microbiome",
}


STUDY_EVIDENCE = {
    "CRA017492": (
        "NMN restored ovarian NAD+, reduced ovarian atrophy and inflammation, increased follicle numbers, and improved hormone secretion and oocyte quality in aged mice.",
        "https://doi.org/10.1002/mco2.727",
    ),
    "GSE128240": (
        "Cyclophosphamide lowered estradiol, increased FSH, reduced ovarian weight, disrupted follicular structure, and downregulated cholesterol-biosynthesis and lipid-metabolism programs.",
        "https://doi.org/10.1016/j.lfs.2019.116999",
    ),
    "GSE194196": (
        "BNC1 deficiency caused premature follicular activation, excessive follicular atresia, oxidative stress, oocyte ferroptosis, and loss of ovarian reserve.",
        "https://doi.org/10.1038/s41467-022-33323-8",
    ),
    "GSE235021": (
        "Whole-grain black rice preserved ovarian reserve and hormonal balance, increased primary and antral follicles, reduced atretic follicles and senescence staining, and partly reversed aging-associated ovarian transcription at 9 months. The central phenotypic validation was specific to black rice.",
        "https://doi.org/10.1007/s11427-025-3046-4",
    ),
    "GSE240547": (
        "High-fat/high-sucrose feeding caused irregular estrous cycles, hormonal disruption, and adverse ovarian transcriptomic and metabolic changes.",
        "https://doi.org/10.1186/s12967-024-04952-y",
    ),
    "GSE249151": (
        "VCD altered sex-hormone levels, follicle numbers, and estrous cycling at days 15 and 30, with stronger day-30 disruption of steroid synthesis, DNA repair, and oocyte-meiosis programs.",
        "https://doi.org/10.7717/peerj.17251",
    ),
    "GSE252791": (
        "Rnf111 haploinsufficiency reduced ovarian reserve and litter size and dysregulated ovarian BMP and gonadal-development pathways.",
        "https://doi.org/10.1007/s11427-024-2606-6",
    ),
    "GSE273724": (
        "HEP14 enhanced mitophagy, reduced oxidative stress, promoted follicular regeneration, restored endocrine function, and shifted the aged ovarian transcriptome toward an adult profile.",
        "https://doi.org/10.1038/s42003-025-08576-w",
    ),
    "GSE279446": (
        "Infertility arose primarily through pituitary Zeb1 derepression, reduced luteinizing-hormone production, and anovulation; ovarian Zeb1 expression was reported as unchanged.",
        "https://doi.org/10.1093/nar/gkaf1357",
    ),
    "GSE279573": (
        "Phytosterol accumulation caused reversible infertility by suppressing folliculogenesis and oocyte quality through ovarian epigenomic changes; a vegetable-oil-free diet rescued fertility.",
        "https://doi.org/10.1038/s42003-024-07233-y",
    ),
    "GSE329024": (
        "CUMS produced age- and duration-dependent HPO-axis responses, with early compensatory changes followed by reduced ovarian size and follicle numbers after longer exposure.",
        "https://doi.org/10.64898/2026.04.24.720585",
    ),
    "GSE98497": (
        "Oocyte-specific Mtor knockout caused infertility; primordial-oocyte knockout produced defective follicular development, progressive oocyte degeneration, loss of granulosa-cell identity, and accumulated oocyte DNA damage.",
        "https://doi.org/10.1073/pnas.1800352115",
    ),
    "PRJNA1076509": (
        "Estropausal-donor fecal microbiota transplantation reduced ovarian inflammation-related and age-upregulated transcription and improved ovarian-health and fertility measures in young recipients.",
        "https://doi.org/10.1038/s43587-026-01069-3",
    ),
    "PRJNA787028": (
        "Ten weeks of high-density housing caused chronic stress and diminished ovarian reserve through granulosa-cell apoptosis, with angiotensinogen/angiotensin-II activation and reduced follicle numbers.",
        "https://doi.org/10.3390/ijms23158614",
    ),
    "PRJNA863443": (
        "The originating study compared natural aging, age-stratified VCD follicle depletion, and Foxl2 haploinsufficiency, identifying shared and model-specific follicle loss, endocrine disruption, and transcriptional remodeling; VCD and Foxl2 models showed distinct hormonal and immune alterations.",
        "https://doi.org/10.1101/2025.06.10.658920",
    ),
    "PRJNA1002222": (
        "CD38 deletion increased ovarian NAD+ levels, preserved follicle reserve and fertility, and partially reversed aging-associated ovarian transcriptional changes at 8 and 12 months.",
        "https://doi.org/10.1038/s43587-023-00532-9",
    ),
}


NO_MONOTONIC_EXPECTATION = {
    "GSE235021_white_brown_vs_polished_6m",
    "GSE235021_white_brown_vs_polished_9m",
    "GSE279446_adult_200H_vs_WT",
    "GSE279446_adult_200M_vs_WT",
    "GSE279446_young_200H_vs_WT",
    "GSE279446_young_200M_vs_WT",
    "GSE329024_CUMS_2M_12W_vs_control",
    "GSE329024_CUMS_2M_8W_vs_control",
    "GSE329024_CUMS_6M_12W_vs_control",
    "GSE329024_CUMS_6M_8W_vs_control",
    "PRJNA863443_Foxl2_het_vs_WT_mid_age",
    "PRJNA863443_Foxl2_het_vs_WT_young",
}


EXPECTED_YOUNGER = {
    "CRA017492_NMN_vs_aged_water",
    "GSE235021_black_vs_polished_6m",
    "GSE235021_black_vs_polished_9m",
    "GSE273724_adult_HEP14_vs_vehicle",
    "GSE273724_aged_HEP14_vs_vehicle",
    "GSE279573_KO_oilfree_vs_KO_chow",
    "PRJNA1076509_estropausal_vs_young_FMT",
    "PRJNA1002222_Cd38KO_vs_WT_8M",
    "PRJNA1002222_Cd38KO_vs_WT_12M",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"No rows for {path.name}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def clean(value: object) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text in {"NA", "NaN", "nan", '""'} else text


def as_int(value: object) -> int:
    text = clean(value)
    if not text:
        return 0
    return int(float(text))


def as_float(value: object) -> float | None:
    text = clean(value)
    if not text:
        return None
    return float(text)


def humanize(value: object) -> str:
    text = clean(value)
    if not text:
        return ""
    return "; ".join(
        re.sub(r"_+", " ", part).strip() for part in text.split(";") if part.strip()
    )


def age_range(low: object, high: object) -> str:
    lo, hi = as_float(low), as_float(high)
    if lo is None or hi is None:
        return "not reported"
    fmt = lambda x: f"{x:.2f}".rstrip("0").rstrip(".")
    return fmt(lo) if math.isclose(lo, hi) else f"{fmt(lo)}–{fmt(hi)}"


def analysis_role(row: dict[str, str]) -> str:
    roles = []
    model_n = as_int(row.get("n_model_ready_samples"))
    projection_n = as_int(row.get("n_external_projection_samples"))
    if model_n:
        split_roles = humanize(row.get("model_ready_split_roles"))
        roles.append(f"frozen benchmark ({split_roles})" if split_roles else "frozen benchmark")
    if projection_n:
        roles.append("external projection")
    if not roles:
        roles.append("metadata inventory")
    return "; ".join(roles)


def build_table1() -> list[dict[str, object]]:
    overview = read_tsv(METHODS / "01_dataset_overview.tsv")
    summary = {
        (row["dataset_folder"], row["modality"]): row
        for row in read_tsv(METHODS / "12_original_dataset_summary.tsv")
    }
    rows = []
    for row in overview:
        dataset = row["dataset_folder"]
        modality = clean(row["modality"])
        source = summary[(dataset, modality)]
        short_ref, title, identifier, status = REFERENCES[dataset]
        rows.append(
            {
                "table_section": (
                    "A. Bulk RNA-seq datasets"
                    if modality == "bulk_rnaseq"
                    else "B. scRNA-derived pseudobulk datasets"
                ),
                "dataset_id": dataset,
                "repository_or_project_accession": clean(row["project_accession"])
                or clean(row["source_accession"])
                or dataset,
                "source_reference": short_ref,
                "source_title": title,
                "doi_or_repository_identifier": identifier,
                "reference_status": status,
                "modality": (
                    "bulk RNA-seq"
                    if modality == "bulk_rnaseq"
                    else "scRNA-derived pseudobulk"
                ),
                "study_design_or_sample_context": humanize(source["sample_contexts"])
                or humanize(row["sample_contexts"]),
                "mouse_strain": humanize(row["strain_canonical"]) or "not reported",
                "age_range_weeks": age_range(
                    source["min_age_week"], source["max_age_week"]
                ),
                "n_curated_samples": as_int(source["n_samples_total"]),
                "n_adult_age_eligible": as_int(
                    source["n_age_range_clock_candidates"]
                ),
                "n_control_like_candidates": as_int(
                    source["n_control_like_clock_candidates"]
                ),
                "n_frozen_benchmark_samples": as_int(
                    row["n_model_ready_samples"]
                ),
                "n_external_projection_samples": as_int(
                    row["n_external_projection_samples"]
                ),
                "analysis_role": analysis_role(row),
                "projection_classes": humanize(row["external_projection_classes"])
                or "not applicable",
                "n_prespecified_tage_contrasts": as_int(row["n_tage_contrasts"]),
                "key_limitation_or_interpretive_note": LIMITATIONS.get(
                    dataset,
                    humanize(row["source_limitations"])
                    or "No dataset-specific limitation beyond the general cross-study harmonization constraints.",
                ),
            }
        )
    return rows


def build_table1_age_bins() -> list[dict[str, object]]:
    """Publish the audited original age-bin summary without redefining bins."""
    rows = read_tsv(METHODS / "15_original_age_bin_summary.tsv")
    return [
        {
            "modality": row["modality"],
            "reproductive_age_bin": row["reproductive_age_bin"],
            "age_range_week": row["age_range_week"],
            "age_range_month": row["age_range_month"],
            "age_range_note": row["age_range_note"],
            "N": as_int(row["N"]),
        }
        for row in rows
    ]


def build_table1_dataset_summary() -> list[dict[str, object]]:
    """Publish the audited original dataset summary without altering its fields."""
    return read_tsv(METHODS / "12_original_dataset_summary.tsv")


TABLE2_COLUMNS = [
    "sample_key",
    "modality",
    "dataset_folder",
    "dataset_name",
    "sample_id",
    "accession",
    "project_accession",
    "sample_title",
    "age_week",
    "age_month_equivalent",
    "age_curation_confidence",
    "bio_age_bin",
    "reproductive_age_bin",
    "strain",
    "strain_canonical",
    "sex",
    "tissue",
    "genotype",
    "treatment",
    "sample_context",
    "developmental_stage",
    "data_format",
    "n_cells_used",
    "n_harmonized_cell_classes",
    "harmonized_cell_classes",
    "age_range_clock_candidate",
    "control_like_clock_candidate",
    "model_building_eligible",
    "projection_reference_control",
    "eligibility_or_exclusion_reason",
    "frozen_split_role",
    "frozen_benchmark_model_ready",
    "external_projection_or_sensitivity",
    "external_projection_boundary",
    "relative_control_limitation",
]


def build_table2() -> list[dict[str, object]]:
    source_rows = read_tsv(METHODS / "02_sample_metadata_publication.tsv")
    rows = []
    for source in source_rows:
        rows.append(
            {
                "sample_key": clean(source["sample_key"]),
                "modality": (
                    "bulk RNA-seq"
                    if source["modality"] == "bulk_rnaseq"
                    else "scRNA-derived pseudobulk"
                ),
                "dataset_folder": clean(source["dataset_folder"]),
                "dataset_name": clean(source["dataset_name"]),
                "sample_id": clean(source["sample_id"]),
                "accession": clean(source["accession"]),
                "project_accession": clean(source["project_accession"]),
                "sample_title": clean(source["sample_title"]),
                "age_week": clean(source["age_week"]),
                "age_month_equivalent": clean(source["age_month_equivalent"]),
                "age_curation_confidence": clean(source["age_curation_confidence"]),
                "bio_age_bin": humanize(source["bio_age_bin"]),
                "reproductive_age_bin": humanize(source["reproductive_age_bin"]),
                "strain": clean(source["strain"]),
                "strain_canonical": humanize(source["strain_canonical"]),
                "sex": clean(source["sex"]),
                "tissue": clean(source["tissue"]),
                "genotype": clean(source["genotype"]),
                "treatment": clean(source["treatment"]),
                "sample_context": humanize(source["sample_context"]),
                "developmental_stage": humanize(source["dev_stage"]),
                "data_format": humanize(source["dataformat"]),
                "n_cells_used": clean(source["n_cells_used"]),
                "n_harmonized_cell_classes": clean(source["n_classes_used"]),
                "harmonized_cell_classes": humanize(
                    source["included_harmonized_cell_classes"]
                ),
                "age_range_clock_candidate": clean(
                    source["age_range_clock_candidate"]
                ),
                "control_like_clock_candidate": clean(
                    source["control_like_clock_candidate"]
                ),
                "model_building_eligible": clean(source["model_building_eligible"]),
                "projection_reference_control": clean(
                    source["projection_reference_control"]
                ),
                "eligibility_or_exclusion_reason": humanize(
                    source["exclusion_reason"]
                ),
                "frozen_split_role": humanize(source["split_role"]),
                "frozen_benchmark_model_ready": clean(
                    source["frozen_v1_model_ready"]
                ),
                "external_projection_or_sensitivity": clean(
                    source["later_external_projection_or_sensitivity"]
                ),
                "external_projection_boundary": humanize(
                    source["projection_benchmark_boundary"]
                ),
                "relative_control_limitation": clean(
                    source["relative_control_limitation"]
                ),
            }
        )
    assert list(rows[0]) == TABLE2_COLUMNS
    return rows


def build_table2_model_development(
    all_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    split_labels = {"train": "training", "validation": "validation", "test": "test"}
    rows = []
    for row in all_rows:
        if row["frozen_benchmark_model_ready"] != "TRUE":
            continue
        split = str(row["frozen_split_role"])
        if split not in split_labels:
            raise RuntimeError(f"Unexpected frozen split role: {split}")
        rows.append(
            {
                "model_development_split": split_labels[split],
                **{
                    **row,
                    "frozen_split_role": split_labels[split],
                },
            }
        )
    split_order = {"training": 0, "validation": 1, "test": 2}
    rows.sort(
        key=lambda row: (
            split_order[row["model_development_split"]],
            str(row["dataset_folder"]),
            str(row["sample_key"]),
        )
    )
    return rows


def build_table2_external_projection(
    all_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    sample_manifest = {
        row["sample_key"]: row
        for row in read_tsv(
            PROJECT
            / "results"
            / "external_projection"
            / "tAge"
            / "manifests"
            / "tage_sample_manifest.tsv"
        )
    }
    external_sample_metadata = {
        row["sample_key"]: row
        for row in read_tsv(METHODS / "04_external_projection_samples.tsv")
    }
    development_datasets = {
        str(row["dataset_folder"])
        for row in all_rows
        if row["frozen_benchmark_model_ready"] == "TRUE"
        and row["frozen_split_role"] in {"train", "validation"}
    }
    heldout_test_datasets = {
        str(row["dataset_folder"])
        for row in all_rows
        if row["frozen_benchmark_model_ready"] == "TRUE"
        and row["frozen_split_role"] == "test"
    }
    rows = []
    for row in all_rows:
        if row["external_projection_or_sensitivity"] != "TRUE":
            continue
        dataset = str(row["dataset_folder"])
        manifest = sample_manifest.get(row["sample_key"])
        if manifest is None:
            external = external_sample_metadata[row["sample_key"]]
            manifest = {
                "projection_class": external["projection_class"],
                "integrated_clock_overlap_status": "fully_external_dataset",
                "developmental_ood": (
                    "TRUE"
                    if "developmental" in external["projection_class"].lower()
                    or "ood" in external["projection_class"].lower()
                    else "FALSE"
                ),
            }
        projection_reference = row["projection_reference_control"] == "TRUE"
        frozen_ready = row["frozen_benchmark_model_ready"] == "TRUE"
        frozen_split = str(row["frozen_split_role"])
        if dataset in development_datasets:
            cohort_group = "cohort overlaps frozen model development"
        elif dataset in heldout_test_datasets and frozen_ready and frozen_split == "test":
            cohort_group = "held-out test cohort also externally projected"
        else:
            cohort_group = "projection-only cohort"
        if projection_reference:
            sample_relationship = (
                "projection-only matched reference control; did not contribute to "
                "model building or model selection"
            )
            overlap_status = "projection only reference control"
        elif frozen_ready and frozen_split == "train":
            sample_relationship = "sample contributed to frozen training"
            overlap_status = humanize(manifest["integrated_clock_overlap_status"])
        elif frozen_ready and frozen_split == "validation":
            sample_relationship = "sample contributed to frozen validation and model selection"
            overlap_status = humanize(manifest["integrated_clock_overlap_status"])
        elif frozen_ready and frozen_split == "test":
            sample_relationship = "sample was held out for final frozen testing"
            overlap_status = humanize(manifest["integrated_clock_overlap_status"])
        else:
            sample_relationship = "projection-only sample; did not contribute to model development"
            overlap_status = humanize(manifest["integrated_clock_overlap_status"])
        rows.append(
            {
                "external_projection_cohort_group": cohort_group,
                "sample_model_development_relationship": sample_relationship,
                "sample_integrated_clock_overlap_status": overlap_status,
                "sample_projection_class": humanize(manifest["projection_class"]),
                "sample_developmental_or_ood": clean(manifest["developmental_ood"]),
                **row,
            }
        )
    rows.sort(
        key=lambda row: (
            str(row["external_projection_cohort_group"]),
            str(row["dataset_folder"]),
            str(row["sample_key"]),
        )
    )
    return rows


def build_table3() -> list[dict[str, object]]:
    source_rows = read_tsv(METHODS / "18_tage_contrast_definitions.tsv")
    rows = []
    for source in source_rows:
        rows.append(
            {
                "dataset_id": clean(source["dataset"]),
                "contrast_id": clean(source["tage_contrast_id"]),
                "contrast_name": clean(source["tage_contrast_name"]),
                "contrast_type": humanize(source["tage_contrast_type"]),
                "reference_group_id": clean(source["tage_reference_group_id"]),
                "n_reference_or_control_samples": as_int(
                    source["n_reference_or_control_samples"]
                ),
                "n_experimental_or_target_samples": as_int(
                    source["n_experimental_or_target_samples"]
                ),
                "n_total_contrast_memberships": as_int(source["n_samples_total"]),
                "age_matching": clean(source["control_matching_age"]),
                "genotype_matching": clean(source["control_matching_genotype"]),
                "treatment_matching": clean(source["control_matching_treatment"]),
                "timing_matching": clean(source["control_matching_timing"]),
                "duration_matching": clean(source["control_matching_duration"]),
                "batch_matching": clean(source["control_matching_batch"]),
                "matching_status": humanize(source["control_matching_status"]),
                "placement_only": clean(source["tage_placement_only"]),
                "developmental_or_ood": (
                    "TRUE"
                    if "developmental" in source["tage_contrast_type"].lower()
                    or "ood" in source["tage_contrast_type"].lower()
                    else "FALSE"
                ),
                "confounding_flag": humanize(source["tage_confounding_flag"]),
                "contrast_warning": clean(source["contrast_warning"]),
                "inclusion_reason": humanize(source["inclusion_reason"]),
                "exclusion_reason": humanize(source["exclusion_reason"]),
                "belongs_to_newly_added_seven": clean(
                    source["external_projection_added7"]
                ),
                "post_lock_prjna1002222_addition": clean(
                    source.get("external_projection_prjna1002222_addition")
                ),
            }
        )
    return rows


def beta_continued_fraction(a: float, b: float, x: float) -> float:
    """Numerical Recipes continued fraction for the incomplete beta function."""
    max_iterations, epsilon, fpmin = 200, 3e-14, 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < fpmin:
        d = fpmin
    d = 1.0 / d
    result = d
    for iteration in range(1, max_iterations + 1):
        m2 = 2 * iteration
        aa = iteration * (b - iteration) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        result *= d * c
        aa = -(a + iteration) * (qab + iteration) * x / (
            (a + m2) * (qap + m2)
        )
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        delta = d * c
        result *= delta
        if abs(delta - 1.0) < epsilon:
            return result
    raise ArithmeticError("Incomplete-beta continued fraction did not converge")


def regularized_incomplete_beta(a: float, b: float, x: float) -> float:
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    front = math.exp(
        math.lgamma(a + b)
        - math.lgamma(a)
        - math.lgamma(b)
        + a * math.log(x)
        + b * math.log1p(-x)
    )
    if x < (a + 1.0) / (a + b + 2.0):
        return front * beta_continued_fraction(a, b, x) / a
    return 1.0 - front * beta_continued_fraction(b, a, 1.0 - x) / b


def welch_t_test(control: list[float], experimental: list[float]) -> float:
    """Two-sided Welch t-test P value, implemented with the standard library."""
    if len(control) < 2 or len(experimental) < 2:
        return math.nan
    mean_c = sum(control) / len(control)
    mean_e = sum(experimental) / len(experimental)
    var_c = sum((value - mean_c) ** 2 for value in control) / (len(control) - 1)
    var_e = sum((value - mean_e) ** 2 for value in experimental) / (
        len(experimental) - 1
    )
    term_c, term_e = var_c / len(control), var_e / len(experimental)
    if math.isclose(term_c + term_e, 0.0):
        return 1.0 if math.isclose(mean_c, mean_e) else 0.0
    t_statistic = (mean_e - mean_c) / math.sqrt(term_c + term_e)
    degrees_freedom = (term_c + term_e) ** 2 / (
        term_c**2 / (len(control) - 1) + term_e**2 / (len(experimental) - 1)
    )
    x = degrees_freedom / (degrees_freedom + t_statistic**2)
    return regularized_incomplete_beta(degrees_freedom / 2.0, 0.5, x)


def student_t_quantile(probability: float, degrees_freedom: float) -> float:
    """Invert the Student-t CDF by bisection using the beta implementation above."""
    if not 0.5 < probability < 1.0 or degrees_freedom <= 0:
        raise ValueError("Student-t quantile requires 0.5 < p < 1 and df > 0")

    def cdf(value: float) -> float:
        x = degrees_freedom / (degrees_freedom + value * value)
        upper_tail_twice = regularized_incomplete_beta(
            degrees_freedom / 2.0, 0.5, x
        )
        return 1.0 - 0.5 * upper_tail_twice

    low, high = 0.0, 1.0
    while cdf(high) < probability:
        high *= 2.0
    for _ in range(100):
        middle = (low + high) / 2.0
        if cdf(middle) < probability:
            low = middle
        else:
            high = middle
    return (low + high) / 2.0


def integrated_effect_stats(
    control: list[float], experimental: list[float]
) -> dict[str, float | int]:
    """Experimental-minus-control effect and Welch 95% confidence interval."""
    n_control, n_experimental = len(control), len(experimental)
    mean_control = sum(control) / n_control
    mean_experimental = sum(experimental) / n_experimental
    var_control = sum((x - mean_control) ** 2 for x in control) / (n_control - 1)
    var_experimental = sum((x - mean_experimental) ** 2 for x in experimental) / (
        n_experimental - 1
    )
    term_control = var_control / n_control
    term_experimental = var_experimental / n_experimental
    standard_error = math.sqrt(term_control + term_experimental)
    degrees_freedom = (term_control + term_experimental) ** 2 / (
        term_control**2 / (n_control - 1)
        + term_experimental**2 / (n_experimental - 1)
    )
    critical = student_t_quantile(0.975, degrees_freedom)
    effect = mean_experimental - mean_control
    return {
        "n_control": n_control,
        "n_experimental": n_experimental,
        "control_mean": mean_control,
        "experimental_mean": mean_experimental,
        "effect_weeks": effect,
        "standard_error": standard_error,
        "confidence_interval_low": effect - critical * standard_error,
        "confidence_interval_high": effect + critical * standard_error,
    }


def bh_adjust(p_values: list[float]) -> list[float]:
    order = sorted(range(len(p_values)), key=p_values.__getitem__)
    adjusted = [math.nan] * len(p_values)
    running = 1.0
    total = len(p_values)
    for rank_index in range(total - 1, -1, -1):
        original_index = order[rank_index]
        rank = rank_index + 1
        running = min(running, p_values[original_index] * total / rank)
        adjusted[original_index] = min(1.0, running)
    return adjusted


def p_text(p_value: float) -> str:
    if p_value < 0.0001:
        return "P<0.0001"
    return f"P={p_value:.4f}"


def clock_result_text(
    elastic_effect: float,
    elastic_p: float,
    pca_effect: float,
    pca_p: float,
) -> str:
    def direction(effect: float) -> str:
        return "younger" if effect < 0 else "older"

    e_sig, p_sig = elastic_p < 0.05, pca_p < 0.05
    same_direction = (elastic_effect < 0) == (pca_effect < 0)
    if e_sig and p_sig and same_direction:
        return (
            f"Both clocks shifted significantly {direction(elastic_effect)} "
            f"(ElasticNet {elastic_effect:+.2f} weeks, {p_text(elastic_p)}; "
            f"PCA/glmnet {pca_effect:+.2f} weeks, {p_text(pca_p)})."
        )
    if e_sig and p_sig:
        return (
            "Both clocks were significant but discordant in direction "
            f"(ElasticNet {elastic_effect:+.2f} weeks, {p_text(elastic_p)}; "
            f"PCA/glmnet {pca_effect:+.2f} weeks, {p_text(pca_p)})."
        )
    if e_sig or p_sig:
        significant_name = "ElasticNet" if e_sig else "PCA/glmnet"
        return (
            f"{significant_name} detected a significant {direction(elastic_effect if e_sig else pca_effect)} shift; "
            f"the other clock was nonsignificant (ElasticNet {elastic_effect:+.2f} weeks, "
            f"{p_text(elastic_p)}; PCA/glmnet {pca_effect:+.2f} weeks, {p_text(pca_p)})."
        )
    relation = (
        f"both shifted {direction(elastic_effect)}"
        if same_direction
        else "the clocks differed in direction"
    )
    return (
        f"Neither clock was significant and {relation} "
        f"(ElasticNet {elastic_effect:+.2f} weeks, {p_text(elastic_p)}; "
        f"PCA/glmnet {pca_effect:+.2f} weeks, {p_text(pca_p)})."
    )


def interpretation_text(
    contrast_id: str,
    elastic_effect: float,
    elastic_p: float,
    pca_effect: float,
    pca_p: float,
) -> tuple[str, str]:
    e_sig, p_sig = elastic_p < 0.05, pca_p < 0.05
    same_direction = (elastic_effect < 0) == (pca_effect < 0)
    if contrast_id in NO_MONOTONIC_EXPECTATION:
        expectation = None
    elif contrast_id in EXPECTED_YOUNGER:
        expectation = -1
    else:
        expectation = 1

    if e_sig and p_sig and same_direction:
        grade = "strong cross-clock support"
    elif (e_sig or p_sig) and same_direction:
        grade = "partial cross-clock support"
    elif not e_sig and not p_sig and same_direction:
        grade = "directional but statistically inconclusive"
    elif e_sig and p_sig:
        grade = "significant cross-clock discordance"
    else:
        grade = "no clock-level replication"

    if expectation is None:
        concordance = (
            "The originating study does not define a single monotonic ovarian-age expectation for this exact contrast, "
            "so the result should be interpreted as a descriptive transcriptomic shift."
        )
    else:
        elastic_match = (elastic_effect * expectation) > 0
        pca_match = (pca_effect * expectation) > 0
        if elastic_match and pca_match:
            concordance = (
                "The direction is phenotype-concordant in both clocks, but statistical strength must be described exactly as reported."
            )
        elif elastic_match or pca_match:
            concordance = (
                "Phenotype concordance is confined to one clock; this is not a two-clock replication."
            )
        else:
            concordance = (
                "The clock direction is not phenotype-concordant as an ovarian-age shift; disease or fertility impairment should not be equated automatically with older transcriptomic age."
            )

    extra = ""
    if contrast_id == "GSE249151_VCD30_vs_vehicle":
        extra = " The day-30 comparison is age-confounded because no age-matched day-30 vehicle was available."
    elif contrast_id.startswith("GSE235021_"):
        extra = " Bulk and scRNA-derived pseudobulk controls were intentionally pooled within dataset and age, making the effect modality/batch-confounded."
    elif contrast_id.startswith("GSE279446_"):
        extra = " This is also a biologically informative negative-control setting because the reported causal lesion was primarily extra-ovarian."
    return grade, concordance + extra


def build_table4() -> list[dict[str, object]]:
    results_dir = PROJECT / "results" / "external_projection"
    membership = read_tsv(METHODS / "19_tage_sample_contrast_membership.tsv")
    predictions = read_tsv(results_dir / "external_projection_predictions_all.tsv")
    effects = read_tsv(
        results_dir
        / "tAge"
        / "cross_clock_concordance"
        / "tables"
        / "cross_clock_contrast_effects.tsv"
    )
    tiers = read_tsv(
        results_dir
        / "tAge"
        / "cross_clock_concordance"
        / "audit"
        / "analysis_tier_membership.tsv"
    )
    table3 = {row["contrast_id"]: row for row in build_table3()}
    primary_ids = {
        contrast_id
        for contrast_id, row in table3.items()
        if row["placement_only"] != "TRUE" and row["developmental_or_ood"] != "TRUE"
    }
    effect_index = {
        (row["tage_contrast_id"], row["clock_name"]): row
        for row in effects
        if row["tage_contrast_id"] in primary_ids
        and row["clock_name"] in {"integrated_elasticnet", "integrated_pca_glmnet"}
    }
    prediction_index = {
        (row["sample_key"], row["model_family"]): float(row["predicted_age_week"])
        for row in predictions
    }
    prediction_rows = {
        (row["sample_key"], row["model_family"]): row for row in predictions
    }
    values: dict[tuple[str, str, str], list[float]] = {}
    for row in membership:
        contrast_id = row["tage_contrast_id"]
        role = row["tage_contrast_role"]
        if contrast_id not in primary_ids or role not in {"control", "experimental"}:
            continue
        for model in ("elasticnet", "pca_glmnet"):
            values.setdefault((contrast_id, model, role), []).append(
                prediction_index[(row["sample_key"], model)]
            )

    # Backward-compatible fallback: derive PRJNA1002222 integrated-clock effect
    # rows only if an older cross-clock table does not yet contain them.
    for contrast_id in sorted(primary_ids):
        for model, clock in (
            ("elasticnet", "integrated_elasticnet"),
            ("pca_glmnet", "integrated_pca_glmnet"),
        ):
            if (contrast_id, clock) in effect_index:
                continue
            if not contrast_id.startswith("PRJNA1002222_"):
                raise RuntimeError(f"Missing integrated effect row for {contrast_id}/{clock}")
            control = values[(contrast_id, model, "control")]
            experimental = values[(contrast_id, model, "experimental")]
            stats = integrated_effect_stats(control, experimental)
            member = next(
                row
                for row in membership
                if row["tage_contrast_id"] == contrast_id
                and row["tage_contrast_role"] == "experimental"
            )
            prediction = prediction_rows[(member["sample_key"], model)]
            effect_index[(contrast_id, clock)] = {
                "dataset": member["dataset"],
                "tage_contrast_id": contrast_id,
                "contrast_name": member["tage_contrast_name"],
                "contrast_type": member["tage_contrast_type"],
                "projection_class": prediction["projection_class"],
                "integrated_clock_overlap_status": "fully_external_dataset",
                **stats,
            }

    working = []
    for contrast_id in sorted(primary_ids):
        elastic = effect_index[(contrast_id, "integrated_elasticnet")]
        pca = effect_index[(contrast_id, "integrated_pca_glmnet")]
        projection_class = elastic["projection_class"]
        if projection_class not in PRIMARY_CLASSES:
            raise RuntimeError(
                f"Unexpected primary projection class for {contrast_id}: {projection_class}"
            )
        elastic_p = welch_t_test(
            values[(contrast_id, "elasticnet", "control")],
            values[(contrast_id, "elasticnet", "experimental")],
        )
        pca_p = welch_t_test(
            values[(contrast_id, "pca_glmnet", "control")],
            values[(contrast_id, "pca_glmnet", "experimental")],
        )
        working.append(
            {
                "contrast_id": contrast_id,
                "elastic": elastic,
                "pca": pca,
                "elastic_p": elastic_p,
                "pca_p": pca_p,
            }
        )

    all_p_values = [
        item[p_name]
        for item in working
        for p_name in ("elastic_p", "pca_p")
    ]
    adjusted = bh_adjust(all_p_values)
    adjusted_iter = iter(adjusted)
    rows = []
    for item in working:
        contrast_id = item["contrast_id"]
        elastic, pca = item["elastic"], item["pca"]
        elastic_p, pca_p = item["elastic_p"], item["pca_p"]
        elastic_q, pca_q = next(adjusted_iter), next(adjusted_iter)
        dataset = elastic["dataset"]
        short_ref, title, identifier, status = REFERENCES[dataset]
        evidence, evidence_url = STUDY_EVIDENCE[dataset]
        grade, interpretation = interpretation_text(
            contrast_id,
            float(elastic["effect_weeks"]),
            elastic_p,
            float(pca["effect_weeks"]),
            pca_p,
        )
        definition = table3[contrast_id]
        rows.append(
            {
                "projection_class": humanize(elastic["projection_class"]),
                "dataset_id": dataset,
                "contrast_id": contrast_id,
                "contrast_name": elastic["contrast_name"],
                "contrast_type": humanize(elastic["contrast_type"]),
                "n_control": as_int(elastic["n_control"]),
                "n_experimental": as_int(elastic["n_experimental"]),
                "integrated_elasticnet_effect_weeks": round(
                    float(elastic["effect_weeks"]), 6
                ),
                "integrated_elasticnet_ci95_low": round(
                    float(elastic["confidence_interval_low"]), 6
                ),
                "integrated_elasticnet_ci95_high": round(
                    float(elastic["confidence_interval_high"]), 6
                ),
                "integrated_elasticnet_welch_p": elastic_p,
                "integrated_elasticnet_bh_q_across_68_tests": elastic_q,
                "integrated_pca_glmnet_effect_weeks": round(
                    float(pca["effect_weeks"]), 6
                ),
                "integrated_pca_glmnet_ci95_low": round(
                    float(pca["confidence_interval_low"]), 6
                ),
                "integrated_pca_glmnet_ci95_high": round(
                    float(pca["confidence_interval_high"]), 6
                ),
                "integrated_pca_glmnet_welch_p": pca_p,
                "integrated_pca_glmnet_bh_q_across_68_tests": pca_q,
                "integrated_clock_result": clock_result_text(
                    float(elastic["effect_weeks"]),
                    elastic_p,
                    float(pca["effect_weeks"]),
                    pca_p,
                ),
                "originating_study_reference": short_ref,
                "originating_study_title": title,
                "originating_study_identifier": identifier,
                "originating_study_reference_status": status,
                "evidence_from_originating_study": evidence,
                "evidence_url": evidence_url,
                "evidence_grade": grade,
                "defensible_interpretation": interpretation,
                "confounding_flag": definition["confounding_flag"],
                "contrast_warning": definition["contrast_warning"],
                "integrated_clock_overlap_status": humanize(
                    elastic["integrated_clock_overlap_status"]
                ),
            }
        )
    if len(rows) != 34 or len({row["contrast_id"] for row in rows}) != 34:
        raise RuntimeError("Supplementary Table 5 must contain 34 unique contrasts.")
    if {row["projection_class"].replace(" ", "_") for row in rows} != PRIMARY_CLASSES:
        raise RuntimeError("Supplementary Table 5 does not contain the five requested classes.")
    return rows


def excel_col(index: int) -> str:
    letters = ""
    while index:
        index, remainder = divmod(index - 1, 26)
        letters = chr(65 + remainder) + letters
    return letters


def xml_text(value: object) -> str:
    text = clean(value)
    return escape(text, {'"': "&quot;"})


def sheet_xml(rows: list[list[object]], freeze_row: int = 1) -> str:
    nrows = len(rows)
    ncols = max(len(row) for row in rows)
    dim = f"A1:{excel_col(ncols)}{nrows}"
    widths = []
    for col in range(ncols):
        values = [clean(row[col]) if col < len(row) else "" for row in rows[:500]]
        width = min(55, max(10, max((len(value) for value in values), default=10) + 2))
        widths.append(
            f'<col min="{col + 1}" max="{col + 1}" width="{width}" customWidth="1"/>'
        )
    row_xml = []
    for r_idx, row in enumerate(rows, start=1):
        cells = []
        for c_idx, value in enumerate(row, start=1):
            ref = f"{excel_col(c_idx)}{r_idx}"
            style = 1 if r_idx == 1 else (2 if r_idx % 2 == 0 else 0)
            if isinstance(value, (int, float)) and not isinstance(value, bool):
                cells.append(f'<c r="{ref}" s="{style}"><v>{value}</v></c>')
            else:
                cells.append(
                    f'<c r="{ref}" s="{style}" t="inlineStr"><is><t xml:space="preserve">'
                    f"{xml_text(value)}</t></is></c>"
                )
        row_xml.append(f'<row r="{r_idx}">{"".join(cells)}</row>')
    pane = (
        f'<pane ySplit="{freeze_row}" topLeftCell="A{freeze_row + 1}" '
        f'activePane="bottomLeft" state="frozen"/>'
    )
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        f'<dimension ref="{dim}"/><sheetViews><sheetView workbookViewId="0">{pane}'
        '</sheetView></sheetViews><sheetFormatPr defaultRowHeight="15"/>'
        f'<cols>{"".join(widths)}</cols><sheetData>{"".join(row_xml)}</sheetData>'
        f'<autoFilter ref="{dim}"/><pageMargins left="0.3" right="0.3" top="0.5" '
        'bottom="0.5" header="0.2" footer="0.2"/></worksheet>'
    )


def write_xlsx(
    path: Path, sheets: list[tuple[str, list[list[object]], int]]
) -> None:
    content_types = [
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">',
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>',
        '<Default Extension="xml" ContentType="application/xml"/>',
        '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>',
        '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>',
        '<Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>',
        '<Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>',
    ]
    for index in range(1, len(sheets) + 1):
        content_types.append(
            f'<Override PartName="/xl/worksheets/sheet{index}.xml" '
            'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        )
    content_types.append("</Types>")

    workbook_sheets = "".join(
        f'<sheet name="{escape(name)}" sheetId="{index}" r:id="rId{index}"/>'
        for index, (name, _, _) in enumerate(sheets, start=1)
    )
    workbook = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        f"<sheets>{workbook_sheets}</sheets></workbook>"
    )
    workbook_rels = [
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">',
    ]
    for index in range(1, len(sheets) + 1):
        workbook_rels.append(
            f'<Relationship Id="rId{index}" '
            'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
            f'Target="worksheets/sheet{index}.xml"/>'
        )
    workbook_rels.append(
        f'<Relationship Id="rId{len(sheets) + 1}" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" '
        'Target="styles.xml"/></Relationships>'
    )
    root_rels = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" '
        'Target="xl/workbook.xml"/><Relationship Id="rId2" '
        'Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" '
        'Target="docProps/core.xml"/><Relationship Id="rId3" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" '
        'Target="docProps/app.xml"/></Relationships>'
    )
    styles = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        '<fonts count="2"><font><sz val="10"/><name val="Arial"/></font>'
        '<font><b/><color rgb="FFFFFFFF"/><sz val="10"/><name val="Arial"/></font></fonts>'
        '<fills count="4"><fill><patternFill patternType="none"/></fill>'
        '<fill><patternFill patternType="gray125"/></fill>'
        '<fill><patternFill patternType="solid"><fgColor rgb="FF1F4E78"/><bgColor indexed="64"/></patternFill></fill>'
        '<fill><patternFill patternType="solid"><fgColor rgb="FFEAF2F8"/><bgColor indexed="64"/></patternFill></fill></fills>'
        '<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>'
        '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
        '<cellXfs count="3"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0" '
        'applyAlignment="1"><alignment vertical="top" wrapText="1"/></xf>'
        '<xf numFmtId="0" fontId="1" fillId="2" borderId="0" xfId="0" applyFill="1" '
        'applyFont="1" applyAlignment="1"><alignment vertical="top" wrapText="1"/></xf>'
        '<xf numFmtId="0" fontId="0" fillId="3" borderId="0" xfId="0" applyFill="1" '
        'applyAlignment="1"><alignment vertical="top" wrapText="1"/></xf></cellXfs>'
        '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>'
        '</styleSheet>'
    )
    today = date.today().isoformat()
    core = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" '
        'xmlns:dc="http://purl.org/dc/elements/1.1/" '
        'xmlns:dcterms="http://purl.org/dc/terms/" '
        'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'
        '<dc:title>Ovarian transcriptomic clock supplementary tables</dc:title>'
        '<dc:creator>Reproducible project table builder</dc:creator>'
        f'<dcterms:created xsi:type="dcterms:W3CDTF">{today}T00:00:00Z</dcterms:created>'
        '</cp:coreProperties>'
    )
    app = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties" '
        'xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">'
        '<Application>Python standard-library XLSX writer</Application></Properties>'
    )
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("[Content_Types].xml", "".join(content_types))
        archive.writestr("_rels/.rels", root_rels)
        archive.writestr("xl/workbook.xml", workbook)
        archive.writestr("xl/_rels/workbook.xml.rels", "".join(workbook_rels))
        archive.writestr("xl/styles.xml", styles)
        archive.writestr("docProps/core.xml", core)
        archive.writestr("docProps/app.xml", app)
        for index, (_, rows, freeze_row) in enumerate(sheets, start=1):
            archive.writestr(
                f"xl/worksheets/sheet{index}.xml", sheet_xml(rows, freeze_row)
            )


def matrix(rows: list[dict[str, object]]) -> list[list[object]]:
    headers = list(rows[0])
    return [headers] + [[row[column] for column in headers] for row in rows]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    all_sample_rows = build_table2()
    table1 = build_table1_dataset_summary()
    table2 = build_table1_age_bins()
    table3a = build_table2_model_development(all_sample_rows)
    table3b = build_table2_external_projection(all_sample_rows)
    table4, table5 = build_table3(), build_table4()

    split_counts = {
        split: sum(row["model_development_split"] == split for row in table3a)
        for split in ("training", "validation", "test")
    }
    projection_datasets = {str(row["dataset_folder"]) for row in table3b}
    projection_cohorts = {
        (str(row["dataset_folder"]), str(row["modality"])) for row in table3b
    }
    development_overlap_cohorts = {
        (str(row["dataset_folder"]), str(row["modality"]))
        for row in table3b
        if row["external_projection_cohort_group"]
        == "cohort overlaps frozen model development"
    }
    heldout_test_projection_cohorts = {
        (str(row["dataset_folder"]), str(row["modality"]))
        for row in table3b
        if row["external_projection_cohort_group"]
        == "held-out test cohort also externally projected"
    }
    projection_only_cohorts = (
        projection_cohorts
        - development_overlap_cohorts
        - heldout_test_projection_cohorts
    )
    frozen_heldout_test_cohorts = {
        (str(row["dataset_folder"]), str(row["modality"]))
        for row in table3a
        if row["model_development_split"] == "test"
    }
    cohort_labels = lambda cohorts: ", ".join(
        f"{dataset} [{modality}]" for dataset, modality in sorted(cohorts)
    )
    if len(table1) != 30:
        raise RuntimeError(f"Expected 30 dataset-summary rows; observed {len(table1)}.")
    if len(table2) != 8 or sum(as_int(row["N"]) for row in table2) != 358:
        raise RuntimeError("Expected eight age-bin rows totaling 358 adult samples.")
    if len(table3a) != 117 or split_counts != {
        "training": 83,
        "validation": 14,
        "test": 20,
    }:
        raise RuntimeError(
            f"Expected 117 samples split 83/14/20; observed {len(table3a)} and {split_counts}."
        )
    if (
        len(table3b) != 300
        or len(projection_datasets) != 20
        or len(projection_cohorts) != 21
    ):
        raise RuntimeError(
            "Expected 300 projection samples across 20 accession-level datasets "
            "and 21 dataset-modality cohorts; observed "
            f"{len(table3b)}, {len(projection_datasets)}, and {len(projection_cohorts)}."
        )
    if (
        len(projection_only_cohorts),
        len(development_overlap_cohorts),
        len(heldout_test_projection_cohorts),
    ) != (17, 4, 0):
        raise RuntimeError(
            "Expected 17 projection-only dataset-modality cohorts, four cohorts "
            "overlapping frozen model development, and no externally projected "
            "held-out test cohort; observed "
            f"{len(projection_only_cohorts)}, {len(development_overlap_cohorts)}, "
            f"and {len(heldout_test_projection_cohorts)}."
        )
    if len(frozen_heldout_test_cohorts) != 2:
        raise RuntimeError(
            "Expected two frozen held-out test cohorts (one bulk RNA-seq and one "
            f"scRNA-derived pseudobulk); observed {len(frozen_heldout_test_cohorts)}."
        )
    projection_reference_rows = [
        row for row in table3b if row["projection_reference_control"] == "TRUE"
    ]
    if len(projection_reference_rows) != 22 or any(
        row["model_building_eligible"] == "TRUE"
        or row["frozen_benchmark_model_ready"] == "TRUE"
        for row in projection_reference_rows
    ):
        raise RuntimeError(
            "Expected 22 projection-only reference controls and none eligible or present "
            "in the frozen model-development cohort."
        )

    t1_path = OUT / "Supplementary_Table_1_original_dataset_summary.tsv"
    t2_path = OUT / "Supplementary_Table_2_original_age_bin_summary.tsv"
    t3a_path = OUT / "Supplementary_Table_3a_model_development_samples.tsv"
    t3b_path = OUT / "Supplementary_Table_3b_external_projection_samples.tsv"
    t4_path = OUT / "Supplementary_Table_4_external_projection_contrasts.tsv"
    t5_path = OUT / "Supplementary_Table_5_primary_integrated_clock_projections.tsv"
    write_tsv(t1_path, table1)
    write_tsv(t2_path, table2)
    write_tsv(t3a_path, table3a)
    write_tsv(t3b_path, table3b)
    write_tsv(t4_path, table4)
    write_tsv(t5_path, table5)

    # Remove superseded reader-facing releases so duplicate table numbers cannot
    # be mistaken for active supplementary tables.
    for stale_name in (
        "Supplementary_Table_1_dataset_characteristics.tsv",
        "Supplementary_Table_2_sample_metadata.tsv",
        "Supplementary_Table_1_original_age_bin_summary.tsv",
        "Supplementary_Table_2a_model_development_samples.tsv",
        "Supplementary_Table_2b_external_projection_samples.tsv",
        "Supplementary_Table_3_external_projection_contrasts.tsv",
        "Supplementary_Table_4_primary_integrated_clock_projections.tsv",
    ):
        stale_path = OUT / stale_name
        if stale_path.exists():
            stale_path.unlink()

    readme_rows = [
        ["Item", "Description"],
        ["Supplementary Table 1", "Dataset summary from methods_tables/12_original_dataset_summary.tsv; 30 dataset-modality rows, including PRJNA1002222."],
        ["Supplementary Table 2", "Reproductive-age-bin definitions and sample counts from methods_tables/15_original_age_bin_summary.tsv; eight modality-by-age-bin rows and 358 adult samples."],
        ["Supplementary Table 3a", "The 117 frozen model-development samples, explicitly labeled as training (83), validation (14), or test (20)."],
        ["Supplementary Table 3b", "The 300 external-projection samples span 20 accession-level datasets and 21 dataset-modality cohorts: 17 projection-only cohorts and four cohorts that overlap frozen model development. The 16 adult CRA003645 samples are retained exclusively in the frozen held-out test and excluded from projection; the six developmental/OOD CRA003645 samples remain descriptive projection-only placements. GSE235021 contributes separate bulk RNA-seq and scRNA-derived pseudobulk cohorts. The 22 matched projection-reference controls are explicitly excluded from model building and model selection."],
        ["Frozen held-out test cohorts", cohort_labels(frozen_heldout_test_cohorts)],
        ["Model-development-overlap projection cohorts", cohort_labels(development_overlap_cohorts)],
        ["Externally projected held-out test cohorts", cohort_labels(heldout_test_projection_cohorts) or "None"],
        ["Projection-only dataset-modality cohorts", cohort_labels(projection_only_cohorts)],
        ["Supplementary Table 4", "Definitions and matching logic for the same 41 validated contrasts used for integrated-clock and tAge projection."],
        ["Supplementary Table 5", "Integrated ElasticNet and PCA/glmnet projections for the 34 primary non-placement, non-developmental/OOD contrasts in intervention, diet, genetic-model, stress-model and microbiome classes, including the 8- and 12-month Cd38-knockout comparisons."],
        ["Supplementary Table 5 statistics", "Effects are experimental-minus-control predicted age in weeks. P values are two-sided Welch t-tests performed separately for each clock; Benjamini-Hochberg q values adjust all 68 contrast-by-clock tests."],
        ["Terminology", "Frozen benchmark denotes samples used in the locked model-development benchmark. External projection denotes post-lock application and does not rewrite the frozen benchmark."],
        ["Interpretation", "Developmental/OOD, placement-only, and confounded contrasts remain explicitly flagged. Clock shifts are descriptive transcriptomic phenotypes, not quantitative biological-age acceleration."],
        ["Generated", date.today().isoformat()],
        ["Source tables", "methods_tables/12_original_dataset_summary.tsv; 15_original_age_bin_summary.tsv; 02_sample_metadata_publication.tsv; 18_tage_contrast_definitions.tsv; external_projection_predictions_all.tsv; tAge sample/contrast manifests; cross-clock effect table and analysis-tier audit"],
    ]
    workbook = OUT / "Ovarian_clock_supplementary_tables.xlsx"
    write_xlsx(
        workbook,
        [
            ("README", readme_rows, 1),
            ("Supplementary Table 1", matrix(table1), 1),
            ("Supplementary Table 2", matrix(table2), 1),
            ("Supplementary Table 3a", matrix(table3a), 1),
            ("Supplementary Table 3b", matrix(table3b), 1),
            ("Supplementary Table 4", matrix(table4), 1),
            ("Supplementary Table 5", matrix(table5), 1),
        ],
    )

    readme = OUT / "README.md"
    readme.write_text(
        "# Publication supplementary tables\n\n"
        "This folder contains five numbered reader-facing supplementary tables, with "
        "Supplementary Table 3 divided into parts 3a and 3b. They are generated from "
        "the audited methods-table sources. The underlying audit tables remain available.\n\n"
        "## Files\n\n"
        "- `Supplementary_Table_1_original_dataset_summary.tsv`: 30 "
        "dataset-modality summary rows reproduced from "
        "`methods_tables/12_original_dataset_summary.tsv`.\n"
        "- `Supplementary_Table_2_original_age_bin_summary.tsv`: the eight original "
        "modality-by-reproductive-age-bin rows (358 adult samples), reproduced from "
        "`methods_tables/15_original_age_bin_summary.tsv`.\n"
        "- `Supplementary_Table_3a_model_development_samples.tsv`: 117 frozen "
        "model-development samples labeled as training (83), validation (14), or test (20).\n"
        "- `Supplementary_Table_3b_external_projection_samples.tsv`: 300 external-projection "
        "samples across 20 accession-level datasets and 21 dataset-modality cohorts: "
        "17 projection-only cohorts and four cohorts overlapping frozen model development. "
        "The 16 adult CRA003645 samples are held-out test samples and are excluded from projection; "
        "the six developmental/OOD CRA003645 samples remain projection-only placements. "
        "GSE235021 contributes separate "
        "bulk RNA-seq and scRNA-derived pseudobulk cohorts. The 22 matched projection-reference "
        "controls are explicitly marked as ineligible for model building or selection; "
        "sample-level overlap and developmental/OOD fields remain explicit.\n"
        "- `Supplementary_Table_4_external_projection_contrasts.tsv`: the same 41 validated "
        "integrated-clock and tAge contrasts, including matching variables, group sizes and OOD/confounding flags.\n"
        "- `Supplementary_Table_5_primary_integrated_clock_projections.tsv`: 34 primary "
        "non-placement, non-developmental/OOD contrasts in intervention, diet, genetic-model, "
        "stress-model and microbiome classes, projected with the integrated ElasticNet and "
        "PCA/glmnet clocks. It includes effect estimates, 95% confidence intervals, two-sided "
        "Welch-test P values, BH-adjusted q values, originating-study evidence and bounded "
        "interpretations.\n"
        "- `Ovarian_clock_supplementary_tables.xlsx`: formatted workbook containing a "
        "README tab and the same five numbered tables (six data sheets because Table 3 "
        "has parts 3a and 3b).\n\n"
        "## Cohort accounting\n\n"
        f"- Frozen held-out test cohorts: {cohort_labels(frozen_heldout_test_cohorts)}.\n"
        f"- Model-development-overlap projection cohorts: {cohort_labels(development_overlap_cohorts)}.\n"
        f"- Externally projected held-out test cohorts: {cohort_labels(heldout_test_projection_cohorts) or 'None'}.\n"
        f"- Projection-only dataset-modality cohorts: {cohort_labels(projection_only_cohorts)}.\n\n"
        "## Interpretation\n\n"
        "Frozen benchmark and post-lock external projection are kept distinct. "
        "Developmental/OOD, placement-only and confounded contrasts are retained and "
        "explicitly labelled. External clock shifts are descriptive transcriptomic "
        "phenotypes rather than quantitative measures of biological-age acceleration. "
        "Supplementary Table 5 reports experimental-minus-control effects; negative values "
        "indicate younger predictions. P values are two-sided Welch t-tests within each "
        "contrast and clock, and q values use Benjamini-Hochberg adjustment across all 68 "
        "contrast-by-clock tests.\n\n"
        "## Rebuild\n\n"
        "Run `python3 scripts/13_build_publication_supplementary_tables.py` from any "
        "working directory.\n",
        encoding="utf-8",
    )
    print(
        f"PASS: table1={len(table1)} dataset rows; table2={len(table2)} age-bin rows; "
        f"table3a={len(table3a)} model samples; table3b={len(table3b)} projection samples; "
        f"table4={len(table4)} contrasts; table5={len(table5)} primary contrasts; "
        f"workbook={workbook}"
    )


if __name__ == "__main__":
    main()
