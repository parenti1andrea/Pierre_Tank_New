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

    def detect_block_type(self, raw_bytes):
        """Return block type string if tag found anywhere in block, else 'PARTICLES'."""
        for tag in (b"RUNH", b"RUNE", b"EVTH", b"EVTE"):
            if tag in raw_bytes:
                return tag.decode("ascii")
        return "PARTICLES"

    def parse_evth(self, floats):
        """Extract primary info from EVTH block."""
        primary_id = int(floats[1])  # particle code
        energy = floats[2]            # GeV
        zenith = floats[10]           # radians
        azimuth = floats[11]          # radians
        return primary_id, energy, zenith, azimuth

    def parse_particles_stream(self, floats, chunk_size=10000):
        """Yield secondary particles in chunks (x,y,z,px,py,pz,t,pid)."""
        n = len(floats) // 7
        particles = floats[: n*7].reshape((n, 7))
        for i in range(0, n, chunk_size):
            yield particles[i:i+chunk_size]

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python corsika_reader.py <CORSIKA binary file>")
        sys.exit(1)

    filename = sys.argv[1]
    reader = CorsikaReader(filename)

    event_idx = -1
    particle_count = 0  # initialize before loop

    for idx, length, raw_bytes, floats in reader.blocks():
        tag = reader.detect_block_type(raw_bytes)

        if tag == "EVTH":
            event_idx += 1
            primary_id, energy, zenith, azimuth = reader.parse_evth(floats)
            print(f"\nEvent {event_idx}:")
            print(f"  Primary: id={primary_id}, energy={energy:.3e} GeV, "
                  f"zenith={zenith:.4f} rad, azimuth={azimuth:.4f} rad")
            particle_count = 0  # reset for new event

        elif tag == "PARTICLES":
            for chunk in reader.parse_particles_stream(floats, chunk_size=10000):
                for p in chunk:
                    x, y, z, px, py, pz, t = p[:7]
                    pid = int(p[6])  # adjust if PID is stored differently
                    print(f"    Particle: pid={pid}, x={x:.3f}, y={y:.3f}, z={z:.3f}, "
                          f"px={px:.3f}, py={py:.3f}, pz={pz:.3f}, t={t:.3f}")
                    particle_count += 1

        elif tag == "EVTE":
            print(f"  Total particles in event: {particle_count}")

