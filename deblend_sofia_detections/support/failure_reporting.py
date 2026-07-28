"""Source-level failure reporting for deblending runs."""

from dataclasses import dataclass
from datetime import datetime


@dataclass(frozen=True)
class SourceFailure:
    """Details of an exception raised while processing one SoFiA source."""

    source_id: str
    cube_name: str
    exception_type: str
    reason: str
    traceback_text: str


def format_failure_report(
    failures,
    requested_count,
    succeeded_count,
    generated_at=None,
):
    """Create a human-readable report for one deblending run."""

    if generated_at is None:
        generated_at = datetime.now().astimezone().isoformat(timespec="seconds")

    lines = [
        "Deblend SoFiA source failure report",
        f"Generated: {generated_at}",
        f"Requested sources: {requested_count}",
        f"Succeeded: {succeeded_count}",
        f"Failed: {len(failures)}",
        "",
    ]

    if not failures:
        lines.append("No source-level failures were recorded.")
        return "\n".join(lines) + "\n"

    for index, failure in enumerate(failures, start=1):
        lines.extend(
            [
                f"Failure {index}",
                f"Source ID: {failure.source_id}",
                f"Cube: {failure.cube_name}",
                f"Exception: {failure.exception_type}",
                f"Reason: {failure.reason}",
                "Traceback:",
                failure.traceback_text.rstrip(),
                "",
            ]
        )

    return "\n".join(lines)


def write_failure_report(
    report_path,
    failures,
    requested_count,
    succeeded_count,
):
    """Overwrite ``report_path`` with the current run's failure report."""

    report = format_failure_report(
        failures,
        requested_count=requested_count,
        succeeded_count=succeeded_count,
    )
    with open(report_path, "w", encoding="utf-8") as report_file:
        report_file.write(report)

    return report_path
