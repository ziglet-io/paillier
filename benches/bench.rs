// Copyright (c) Jeffrey Hohenstein <jeffrey.hohenstein@gmail.com>
//
// All rights reserved.

use criterion::{criterion_group, criterion_main, Criterion};
use crypto_bigint::{rand_core::OsRng, RandomMod, Uint};
use paillier::Paillier;

fn bench(criterion: &mut Criterion) {
    let mut group = criterion.benchmark_group("paillier");

    const T: u32 = 10;
    const SINGLE: usize = 16;
    const DOUBLE: usize = 32;
    const QUAD: usize = 64;
    const OCT: usize = 128;

    type P<'a> = Paillier<'a, SINGLE, DOUBLE, QUAD, OCT, OsRng>;

    let mut rng = OsRng;
    let mut rng2 = OsRng;
    let kp = P::generate_key(T, &mut rng);
    let mut paillier = P::new(&kp, &mut rng);
    let message = Uint::<DOUBLE>::random_mod(&mut rng2, &kp.public_key.n.to_nz().unwrap());
    let p = kp
        .secret_key
        .p
        .concat(&Uint::<SINGLE>::ZERO)
        .concat(&Uint::<DOUBLE>::ZERO);
    let q = kp
        .secret_key
        .p
        .concat(&Uint::<SINGLE>::ZERO)
        .concat(&Uint::<DOUBLE>::ZERO);
    let n = p * q;
    let nsquared = n.wrapping_square();

    group.bench_function("random_znsquared 2048bit", |b| {
        b.iter(|| {
            let _r = P::random_znsquared(kp.secret_key.p, kp.secret_key.q, n, nsquared, &mut rng2);
        });
    });

    group.bench_function("encrypt 2048bit message", |b| {
        b.iter(|| {
            let _c = paillier.encrypt(&message);
        });
    });

    group.bench_function("decrypt 2048bit message", |b| {
        let c = paillier.encrypt(&message);
        b.iter(|| paillier.decrypt(&c));
    });
}

criterion_group!(benches, bench);
criterion_main!(benches);
