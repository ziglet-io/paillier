// Copyright (c) Jeffrey Hohenstein <jeffrey.hohenstein@gmail.com>
//
// All rights reserved.

use crypto_bigint::{rand_core::OsRng, Int, Random, RandomMod, Uint};

use crate::Paillier;

const T: u32 = 10;
const SINGLE: usize = 16;
const DOUBLE: usize = 32;
const QUAD: usize = 64;
const OCT: usize = 128;
type P<'a> = Paillier<'a, SINGLE, DOUBLE, QUAD, OCT, OsRng>;

#[test]
fn generate_probable_prime() {
    let mut rng = OsRng;
    let candidate = P::generate_probable_prime(T, &mut rng);
    println!("candidate {}", candidate);
}

#[test]
fn generate_key() {
    let mut rng = OsRng;
    let kp = P::generate_key(T, &mut rng);
    println!("p {}", kp.secret_key.p);
    println!("q {}", kp.secret_key.q);
    println!("n {}", kp.public_key.n);
}

#[test]
fn encrypt() {
    let mut rng = OsRng;
    let mut rng3 = OsRng;
    let kp = P::generate_key(T, &mut rng);
    let mut paillier = P::new(&kp, &mut rng);
    let message = Uint::<DOUBLE>::random_mod(&mut rng3, &kp.public_key.n.to_nz().unwrap());
    let _c = paillier.encrypt(&message);
}

#[test]
fn decrypt() {
    let mut rng = OsRng;
    let kp = P::generate_key(T, &mut rng);
    let message = Uint::<QUAD>::random(&mut rng);
    let paillier = P::new(&kp, &mut rng);
    paillier.decrypt(&message);
}

#[test]
fn encrypt_decrypt() {
    let mut rng = OsRng;
    let kp = P::generate_key(T, &mut rng);
    let message = Uint::<DOUBLE>::random_mod(&mut rng, &kp.public_key.n.to_nz().unwrap());
    println!("message {}", message);
    let mut paillier = P::new(&kp, &mut rng);
    let c = paillier.encrypt(&message);
    println!("c {}", c);
    let m = paillier.decrypt(&c);
    println!("m {}", m);
    assert_eq!(message, m);
}

#[test]
fn gcd() {
    let a = Int::<4>::from_i64(240);
    let b = Int::<4>::from_i64(46);
    let (gcd, _, _) = P::scalar_gcd(a, b);
    assert!(gcd == Int::<4>::from_i64(2));
}

#[test]
fn lambda() {
    let a = Uint::<SINGLE>::from_u64(23);
    let b = Uint::<SINGLE>::from_u64(17);
    let lambda = P::lambda(a, b);
    assert_eq!(lambda.as_words()[0], 176);
}

#[test]
fn random_zn_is_invertible_mod_zn() {
    let mut rng = OsRng;
    let kp = P::generate_key(T, &mut rng);
    let n = kp.public_key.n.concat(&Uint::<DOUBLE>::ZERO);

    for _ in 0..100 {
        let random_n = P::random_zn(kp.secret_key.p, kp.secret_key.q, &mut rng);
        let _ = random_n.inv_mod(&n).expect("inverse");
    }
}

#[test]
fn random_znsquared_is_invertible_mod_znsquared() {
    let mut rng = OsRng;
    let kp = P::generate_key(T, &mut rng);
    let n = kp.public_key.n.concat(&Uint::<DOUBLE>::ZERO);
    let nsquared = n.wrapping_square();

    for _ in 0..100 {
        let random_nsquared =
            P::random_znsquared(kp.secret_key.p, kp.secret_key.q, n, nsquared, &mut rng);
        let _ = random_nsquared
            .inv_mod(&n.wrapping_square())
            .expect("inverse");
    }
}

#[test]
fn totient() {
    let mut rng = OsRng;
    let kp = P::generate_key(T, &mut rng);
    let paillier = P::new(&kp, &mut rng);
    let p = kp
        .secret_key
        .p
        .concat(&Uint::<SINGLE>::ZERO)
        .concat(&Uint::<DOUBLE>::ZERO);
    let q = kp
        .secret_key
        .q
        .concat(&Uint::<SINGLE>::ZERO)
        .concat(&Uint::<DOUBLE>::ZERO);
    let t2 = (p - Uint::ONE) * (q - Uint::ONE);
    assert_eq!(t2, paillier.totient);
}
