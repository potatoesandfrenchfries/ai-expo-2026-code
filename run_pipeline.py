"""Knotworking Drug Discovery Pipeline — Phase 1 entry point.

Chains Layers 2–10 for an end-to-end run from patient record to ranked
drug candidates.

Quick start (no internet or API keys required):
    python run_pipeline.py --demo

Full run (requires network access + ANTHROPIC_API_KEY for Rocq proof):
    ANTHROPIC_API_KEY=sk-ant-... python run_pipeline.py

Options:
    --demo          Use stub data for L3–L5 (no network / API key required).
    --candidates N  Number of molecules to generate (default: 50).
    --output PATH   JSON output file path (default: pipeline_output.json).

The pipeline prints a live progress banner as each layer completes and
produces a unified summary report at the end.  All results are also
serialised to a JSON file for downstream consumption or audit trails.
"""

import argparse

from src.pipeline.pipeline_runner import DrugDiscoveryPipeline, demo_patient


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Knotworking AI Drug Discovery Pipeline — Phase 1",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Run with stub L3–L5 data (no network or API key required).",
    )
    parser.add_argument(
        "--candidates",
        type=int,
        default=50,
        metavar="N",
        help="Number of molecular candidates to generate (default: 50).",
    )
    parser.add_argument(
        "--output",
        default="pipeline_output.json",
        metavar="PATH",
        help="Path for the JSON results file (default: pipeline_output.json).",
    )
    args = parser.parse_args()

    patient  = demo_patient()
    pipeline = DrugDiscoveryPipeline(
        n_candidates=args.candidates,
        demo_mode=args.demo,
    )

    result = pipeline.run(patient)
    pipeline.print_full_report(result)

    output_path = args.output
    if output_path == "pipeline_output.json":
        output_path = f"pipeline_output_{patient.patient_id}.json"
    pipeline.save_json(result, output_path)


if __name__ == "__main__":
    main()
