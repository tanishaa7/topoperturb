"""Download verified prepared tables only; never fetch the large raw count files."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
from urllib.request import urlopen

REPO = "huggingbahl21/causalbench"
REVISION = "068c7902f1740b8d653033c3decfed2c7156fa3d"
SHA256 = {
    "k562": "455024efc9f40b877f7f9c1d55f9eaaf957c7995c6289ffbcd65251fc856fee5",
    "rpe1": "6ddd100eff9dbcab65ff94b662a87e9eb50de6401b5242c0f0108692c92cc408",
}
ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fetch(name, destination):
    destination.mkdir(parents=True, exist_ok=True)
    filename = f"causalbench_{name}.tsv.gz"
    path = destination / filename
    url = f"https://huggingface.co/datasets/{REPO}/resolve/{REVISION}/data/prepared/{filename}"
    if not path.exists():
        partial = path.with_suffix(path.suffix + ".part")
        with urlopen(url, timeout=120) as response, partial.open("wb") as output:
            shutil.copyfileobj(response, output)
        if sha256(partial) != SHA256[name]:
            raise ValueError(f"Checksum mismatch: {partial}; final file was not created")
        partial.replace(path)
    if sha256(path) != SHA256[name]:
        raise ValueError(f"Existing file does not match the pinned source: {path}")
    record = {"repository": REPO, "revision": REVISION, "url": url,
              "sha256": SHA256[name], "bytes": path.stat().st_size}
    path.with_name(filename + ".source.json").write_text(json.dumps(record, indent=2) + "\n")
    print(path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cell_line", choices=["k562", "rpe1", "all"])
    parser.add_argument("--output", type=Path, default=ROOT / "data/cache/causalbench")
    args = parser.parse_args()
    for name in SHA256 if args.cell_line == "all" else [args.cell_line]:
        fetch(name, args.output)
