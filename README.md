# CSI‑FiSh – Isogeny‑Based Signatures in Rust

> **⚠️  Experimental code. Use at your own risk.**

CSI‑FiSh (`csifish` crate) is a pure‑Rust implementation of the CSI-FiSh: Efficient Isogeny based Signatures through Class Group Computations

The design combines **Supersingular Isogeny** techniques with class‑group computations to obtain post‑quantum signatures with very small public keys.

This repository provides:

* Blinded, constant‑time arithmetic for the quadratic‑imaginary class group underlying CSI‑FiSh.
* Highly optimised modular inversion written in hand‑tuned assembly (`x86_64`, `aarch64`) and automatically selected at build time.
* Trait‑based API for key generation, signing and verification.
* Optional **Rayon** parallelisation (`parallel` Cargo feature).

---

## Minimum requirements

|          | Version / tool                                          |
| -------- | ------------------------------------------------------- |
| **Rust** | nightly `1.78` or newer (uses `generic_const_exprs`)    |
| **GMP**  | `6.2+` development headers (via `rug` / `gmp-mpfr-sys`) |
| **OS**   | Linux, macOS, Windows (x86\_64 / AArch64 tested)        |

```bash
# Debian / Ubuntu
sudo apt install libgmp-dev build-essential pkg-config

# Fedora
sudo dnf install gmp-devel
```

---

## Building

```bash
git clone https://github.com/your‑org/csifish.git
cd csifish
rustup toolchain install nightly
rustup override set nightly
cargo build --release        # add --features parallel for multithreading
```

Running the small self‑tests:

```bash
cargo test --release
```

---

## Quick start

```rust
use csifish::csifish::signature::SigningKey;

const CURVES: u32 = 256;   // security parameter
const ROUNDS: u32 = 7;     // number of Fiat–Shamir rounds
const HASHES: u32 = 11;    // leaves per Merkle proof

fn main() {
    // 1.  Generate a signing key (= secret key)
    let sk = SigningKey::<CURVES, ROUNDS, HASHES>::generate();
    let vk = sk.verifying_key();           // derive verifying (public) key

    // 2.  Sign arbitrary messages
    let msg = b"🐟 post‑quantum ahoy!";
    let sig = sk
        .try_sign(msg)                    // Signer trait provides try_sign()
        .expect("failed to sign");

    // 3.  Verify signatures
    vk.verify(msg, &sig).expect("invalid signature");
}
```

See [`examples/`](examples) for complete runnable programs.

---

## Crate features

| Feature           | Default | Description                                                      |
| ----------------- | ------- | ---------------------------------------------------------------- |
| `parallel`        | ❌       | Enable Rayon‑backed `ParallelIterator` implementations           |
| `use-system-libs` | ❌       | Link against the system‑installed GMP instead of the bundled one |
---

## Security Notice

Although care was taken to follow constant‑time coding practices, **this code has not been formally audited**.
Do **not** deploy in production systems or handle sensitive key material until a thorough cryptographic review has been completed.

---

## License

This project is distributed under the terms of the **MIT** license.
See [LICENSE](LICENSE) for details.
