import struct
import numpy as np
import sys

class CorsikaReader:
    def __init__(self, filename):
        self.filename = filename

    def blocks(self):
        """Generator that yields (block_index, block_length, raw_bytes, np.float32 array)."""
        with open(self.filename, "rb") as f:
            idx = 0
            while True:
                length_bytes = f.read(4)
                if not length_bytes:
                    break  # EOF
                length = struct.unpack("i", length_bytes)[0]
                data = f.read(length)
                if len(data) < length:
                    break
                end_length_bytes = f.read(4)
                if not end_length_bytes:
                    break
                end_length = struct.unpack("i", end_length_bytes)[0]
                if length != end_length:
                    raise ValueError(f"Mismatched record markers in block {idx}: {length} vs {end_length}")
                floats = np.frombuffer(data, dtype=np.float32)
                yield idx, length, data, floats
                idx += 1

    def parse_evth(self, floats):
        """Extract primary info from EVTH block."""
        primary_id = int(floats[1])  # particle code
        energy = floats[2]  # GeV
        zenith = floats[10]  # radians
        azimuth = floats[11]  # radians
        return primary_id, energy, zenith, azimuth

    def parse_particles(self, floats):
        """Parse secondary particles from particle data block."""
        n = len(floats) // 7
        particles = floats[: n * 7].reshape((n, 7))
        return particles

    def detect_block_type(self, raw_bytes):
        """Return block type string if tag found anywhere in block, else 'PARTICLES'."""
        for tag in (b"RUNH", b"RUNE", b"EVTH", b"EVTE"):
            if tag in raw_bytes:
                return tag.decode("ascii")
        return "PARTICLES"

def scan_for_tags(filename, tags=(b"RUNH", b"RUNE", b"EVTH", b"EVTE"), max_hits=50):
    with open(filename, "rb") as f:
        data = f.read()
    for tag in tags:
        idx = data.find(tag)
        hit_count = 0
        while idx != -1 and hit_count < max_hits:
            print(f"Found {tag.decode()} at byte offset {idx}")
            hit_count += 1
            idx = data.find(tag, idx + 1)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python corsika_reader.py <CORSIKA binary file> [--list-tags] [--events] [--scan-raw]")
        sys.exit(1)

    filename = sys.argv[1]
    list_tags = "--list-tags" in sys.argv
    summarize_events = "--events" in sys.argv
    scan_raw = "--scan-raw" in sys.argv

    reader = CorsikaReader(filename)

    if list_tags:
        print("Scanning file for block tags:")
        for idx, length, raw_bytes, floats in reader.blocks():
            tag = reader.detect_block_type(raw_bytes)
            if tag != "PARTICLES":
                print(f"Block {idx}: tag={tag}")
        sys.exit(0)

    if summarize_events:
        event_idx = -1
        particle_blocks = 0
        particle_count = 0
        for idx, length, raw_bytes, floats in reader.blocks():
            tag = reader.detect_block_type(raw_bytes)

            if tag == "EVTH":
                event_idx += 1
                primary_id, energy, zenith, azimuth = reader.parse_evth(floats)
                print(f"Event {event_idx}:")
                print(f"  EVTH: primary={primary_id}, energy={energy:.3e} GeV, zenith={zenith:.4f} rad, azimuth={azimuth:.4f} rad")
                particle_blocks = 0
                particle_count = 0

            elif tag == "PARTICLES":
                particles = reader.parse_particles(floats)
                particle_blocks += 1
                particle_count += len(particles)

            elif tag == "EVTE":
                print(f"  {particle_blocks} particle blocks containing {particle_count} particles")
                print("  EVTE")

        sys.exit(0)

    if scan_raw:
        print("Raw scan for tags:")
        scan_for_tags(filename)
        sys.exit(0)

    # Default: just show first few blocks with their type
    for idx, length, raw_bytes, floats in reader.blocks():
        tag = reader.detect_block_type(raw_bytes)
        print(f"Block {idx}: {tag}, length={length}")
        if idx >= 20:
            break

