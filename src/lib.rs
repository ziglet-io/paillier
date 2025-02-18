// Copyright (c) Jeffrey Hohenstein <jeffrey.hohenstein@gmail.com>
//
// All rights reserved.

use crypto_bigint::{
    modular::{MontyForm, MontyParams},
    rand_core::RngCore,
    Concat, Int, InvMod, NonZero, RandomMod, Split, Uint, Zero,
};
use crypto_primality::miller_rabin;

/// Secret Key ($p$,$q$) where $p,q \in \mathbb{Z}_n$ are primes
///
/// # Parameters
/// * `LIMBS_SINGLE`: size of $p,q$ in [Limb]
pub struct SecretKey<const LIMBS_SINGLE: usize> {
    pub p: Uint<LIMBS_SINGLE>,
    pub q: Uint<LIMBS_SINGLE>,
}

/// Public key $n \eq $pq$
/// # Parameters
/// * `LIMBS_DOUBLE`: size of $n \eq pq$ in [Limb]
pub struct PublicKey<const LIMBS_DOUBLE: usize> {
    pub n: Uint<LIMBS_DOUBLE>,
}

/// KeyPair
///
/// # Parameters
/// * `LIMBS_SINGLE`: size of $p,q$ in [Limb]
/// * `LIMBS_DOUBLE`: size of $n \eq pq$ and message in [Limb]
/// * `LIMBS_QUAD`: size of ciphertext in [Limb]
pub struct KeyPair<const LIMBS_SINGLE: usize, const LIMBS_DOUBLE: usize, const LIMBS_QUAD: usize> {
    pub public_key: PublicKey<LIMBS_DOUBLE>,
    pub secret_key: SecretKey<LIMBS_SINGLE>,
}

/// The Paillier Encryption system for a single message of fixed size.
///
/// Implemented using [crypto_bigint] enabling variable sizing of keys, messages, and ciphertexts.
///
/// # Parameters
/// * `LIMBS_SINGLE`: size of $p,q$ in [Limb]
/// * `LIMBS_DOUBLE`: size of $n \eq pq$ and message in [Limb]
/// * `LIMBS_QUAD`: size of ciphertext in [Limb]
/// * `LIMBS_OCT`: twice the size of `LIMBS_QUAD`. Placeholder
/// * `R`: type of random number generator
///
/// # Example
///
/// ```rust
/// use crypto_bigint::{rand_core::OsRng, Int, Random, RandomMod, Uint};
/// use paillier::Paillier;
///
/// const T: u32 = 10;
/// const SINGLE: usize = 16;
/// const DOUBLE: usize = 32;
/// const QUAD: usize = 64;
/// const OCT: usize = 128;
/// type P<'a> = Paillier<'a, SINGLE, DOUBLE, QUAD, OCT, OsRng>;
///
/// let mut rng = OsRng;
/// let mut rng2 = OsRng;
/// let kp = P::generate_key(T, &mut rng);
/// let mut paillier = P::new(&kp, &mut rng2);
///
/// let message = Uint::<DOUBLE>::random_mod(&mut rng, &kp.public_key.n.to_nz().unwrap());
/// let c = paillier.encrypt(&message);
/// let m = paillier.decrypt(&c);
/// ```
///
/// # References
/// * [1](https://link.springer.com/content/pdf/10.1007%2F3-540-48910-X_16.pdf) Paillier, Pascal. “Public-Key Cryptosystems Based on Composite Degree Residuosity Classes.” In Advances in Cryptology — EUROCRYPT ’99, edited by Jacques Stern, 223–38. Berlin, Heidelberg: Springer, 1999. https://doi.org/10.1007/3-540-48910-X_16.
/// * [2](https://www.amazon.com/Introduction-Cryptography-Chapman-Network-Security/dp/0815354363/ref=sr_1_1?sr=8-1) Katz, Jonathan, and Yehuda Lindell. Introduction to Modern Cryptography. 3rd ed. New York: Chapman and Hall/CRC, 2020. https://doi.org/10.1201/9781351133036.
/// * [3](https://github.com/RustCrypto/crypto-bigint) [RustCrypto](https://github.com/RustCrypto/crypto-bigint)
pub struct Paillier<
    'a,
    const LIMBS_SINGLE: usize,
    const LIMBS_DOUBLE: usize,
    const LIMBS_QUAD: usize,
    const LIMBS_OCT: usize,
    R,
> where
    R: RngCore,
{
    /// KeyPair
    key_pair: &'a KeyPair<LIMBS_SINGLE, LIMBS_DOUBLE, LIMBS_QUAD>,
    /// Random number generator for $r$
    rng: &'a mut R,
    /// Widened $n$
    n: Uint<LIMBS_QUAD>,
    /// Precomputed $n^2$
    nsquared: Uint<LIMBS_QUAD>,
    /// Precomputed $\phi(n) \eq (p-1)(q-1)$
    totient: Uint<LIMBS_QUAD>,
}

impl<
        'a,
        const LIMBS_SINGLE: usize,
        const LIMBS_DOUBLE: usize,
        const LIMBS_QUAD: usize,
        const LIMBS_OCT: usize,
        R,
    > Paillier<'a, LIMBS_SINGLE, LIMBS_DOUBLE, LIMBS_QUAD, LIMBS_OCT, R>
where
    Uint<LIMBS_SINGLE>: Concat<Output = Uint<LIMBS_DOUBLE>> + InvMod<Output = Uint<LIMBS_SINGLE>>,
    Uint<LIMBS_DOUBLE>: Concat<Output = Uint<LIMBS_QUAD>> + InvMod<Output = Uint<LIMBS_DOUBLE>>,
    Uint<LIMBS_QUAD>: Concat<Output = Uint<LIMBS_OCT>> + InvMod<Output = Uint<LIMBS_QUAD>>,
    Uint<LIMBS_OCT>: Split<Output = Uint<LIMBS_QUAD>>,
    Uint<LIMBS_QUAD>: Split<Output = Uint<LIMBS_DOUBLE>>,
    Uint<LIMBS_DOUBLE>: Split<Output = Uint<LIMBS_SINGLE>>,
    R: RngCore,
{
    /// Given a reference to an existing [KeyPair], create a new instance of the encryptor/decryptor
    pub fn new(
        key_pair: &'a KeyPair<LIMBS_SINGLE, LIMBS_DOUBLE, LIMBS_QUAD>,
        rng: &'a mut R,
    ) -> Self {
        let p = key_pair
            .secret_key
            .p
            .concat(&Uint::<LIMBS_SINGLE>::ZERO)
            .concat(&Uint::<LIMBS_DOUBLE>::ZERO);
        let q = key_pair
            .secret_key
            .q
            .concat(&Uint::<LIMBS_SINGLE>::ZERO)
            .concat(&Uint::<LIMBS_DOUBLE>::ZERO);

        let totient = (p - Uint::ONE) * (q - Uint::ONE);

        let n = key_pair.public_key.n.concat(&Uint::<LIMBS_DOUBLE>::ZERO);

        let n_squared = n.wrapping_square();

        Self {
            key_pair,
            rng,
            n,
            nsquared: n_squared,
            totient,
        }
    }

    /// Given a message of fixed length encrypt it.
    ///
    /// # Parameters
    /// * `message`: $m \in \mathbb{Z}_{n^2}^\star$
    /// * `LIMBS_DOUBLE`: the size in [Limb] of $m$
    /// * `LIMBS_QUAD`: the size of $n^2$ in [Limb].
    pub fn encrypt(&mut self, message: &Uint<LIMBS_DOUBLE>) -> Uint<LIMBS_QUAD> {
        let r = Self::random_znsquared(
            self.key_pair.secret_key.p,
            self.key_pair.secret_key.q,
            self.n,
            self.nsquared,
            &mut self.rng,
        );

        let r = MontyForm::new(&r, MontyParams::new(self.nsquared.to_odd().unwrap()));

        #[allow(non_snake_case)]
        let N = MontyForm::new(&self.n, MontyParams::new(self.nsquared.to_odd().unwrap()));

        let one = MontyForm::new(
            &Uint::<LIMBS_QUAD>::ONE,
            MontyParams::new(self.nsquared.to_odd().unwrap()),
        );

        let ciphertext = (one + N).pow(message) * r.pow(&self.n);

        let c = ciphertext.retrieve();

        c
    }

    /// Given a ciphertext, decrypt it.
    ///
    /// # Parameters
    /// `ciphertext`: $c \in \mathbb{Z}_{n^2}^\star$
    pub fn decrypt(&self, ciphertext: &Uint<LIMBS_QUAD>) -> Uint<LIMBS_DOUBLE> {
        let c = MontyForm::new(
            &ciphertext,
            MontyParams::new(self.nsquared.to_odd().unwrap()),
        );
        let c_prime = c.pow(&self.totient);

        let m_temp = (c_prime.retrieve() - Uint::ONE) / self.n;
        let m_prime = MontyForm::new(&m_temp, MontyParams::new(self.n.to_odd().unwrap()));

        let totient_inv = MontyForm::new(
            &(self.totient.inv_mod(&self.n).unwrap()),
            MontyParams::new(self.n.to_odd().unwrap()),
        );

        let m = m_prime * totient_inv;

        m.retrieve().split().0
    }

    /// Given two non-zero integers $a$, $b$ compute the greatest common divisor
    /// and also return Bezout's coefficients $s$ and $t$.
    pub fn scalar_gcd<const LIMBS: usize>(
        a: Int<LIMBS>,
        b: Int<LIMBS>,
    ) -> (Int<LIMBS>, Int<LIMBS>, Int<LIMBS>) {
        assert!(a.is_zero().unwrap_u8() == 0u8);
        assert!(b.is_zero().unwrap_u8() == 0u8);

        let mut q: Int<LIMBS>;

        let mut r_minus_1 = if a > b { a } else { b };
        let mut r = if a > b { b } else { a };

        let mut s_minus_1 = Int::<LIMBS>::ONE;
        let mut s = Int::<LIMBS>::ZERO;

        let mut t_minus_1 = Int::<LIMBS>::ZERO;
        let mut t = Int::<LIMBS>::ONE;

        loop {
            let (qx, _rx) = r_minus_1.checked_div_rem(&r.to_nz().unwrap());
            q = qx.expect("division error");

            (r_minus_1, r) = (r, r_minus_1 - q * r);
            (s_minus_1, s) = (s, s_minus_1 - q * s);
            (t_minus_1, t) = (t, t_minus_1 - q * t);

            if r.is_zero().unwrap_u8() == 1u8 {
                break;
            }
        }

        (q, s, t)
    }

    /// Generate a random element of $\mathbb{Z}^\star_{n^2}$
    ///
    /// * $f(a,b) = (1 + n)^ab^n \mod n^2$
    /// * a is an element of $Z_n$
    /// * b is an element of $Z_{n^2}^\star$
    pub fn random_znsquared(
        p: Uint<LIMBS_SINGLE>,
        q: Uint<LIMBS_SINGLE>,
        n: Uint<LIMBS_QUAD>,
        nsquared: Uint<LIMBS_QUAD>,
        mut rng: &mut R,
    ) -> Uint<LIMBS_QUAD> {
        let a = Uint::<LIMBS_QUAD>::random_mod(&mut rng, &n.to_nz().unwrap());
        let b = Self::random_zn(p, q, &mut rng);
        #[allow(non_snake_case)]
        let B = MontyForm::new(&b, MontyParams::new(nsquared.to_odd().unwrap()));
        let one = MontyForm::new(&Uint::ONE, MontyParams::new(nsquared.to_odd().unwrap()));
        #[allow(non_snake_case)]
        let N = MontyForm::new(&n, MontyParams::new(nsquared.to_odd().unwrap()));

        let x = (N + one).pow(&a) * B.pow(&n);

        x.retrieve()
    }

    /// Given two integers $a$ and $b$, compute the least common multiple as
    /// $\frac{\mid{ab}\mid}{gcd(a,b)}$
    pub fn scalar_lcm<const LIMBS: usize>(a: Int<LIMBS>, b: Int<LIMBS>) -> Uint<LIMBS> {
        let (gcd, _, _) = Self::scalar_gcd(a, b);
        let lcm = (a * b).checked_div(&gcd).expect("division error").abs();
        lcm
    }

    /// Compute Carmichael's Number $\lambda(n) = lcm((p-1)(q-1))$ where $n$ is the product of two primes $p$ and $q$
    pub fn lambda<const LIMBS: usize>(p: Uint<LIMBS>, q: Uint<LIMBS>) -> Uint<LIMBS> {
        let p = p.as_int();
        let q = q.as_int();
        let p_minus_1 = p - Int::ONE;
        let q_minus_1 = q - Int::ONE;
        Self::scalar_lcm::<LIMBS>(p_minus_1, q_minus_1)
    }

    /// Attempt to generate a prime number between $0..$ [Uint::<LIMBS_SINGLE>::BITS] - 1
    ///
    /// Uses the [Miller Rabin](miller_rabin::is_composite) test to check that the generated
    /// value is prime with high probability.
    pub fn generate_probable_prime(t: u32, mut rng: &mut R) -> Uint<LIMBS_SINGLE> {
        miller_rabin::generate_probable_prime(Uint::<LIMBS_SINGLE>::BITS, t, &mut rng)
    }

    /// Select an element $x \in \mathbb{Z}_{n}^{*}$ where $n = pq$
    ///
    /// Both elements mod p or q must be coprime to both p and q
    pub fn random_zn(
        p: Uint<LIMBS_SINGLE>,
        q: Uint<LIMBS_SINGLE>,
        mut rng: &mut R,
    ) -> Uint<LIMBS_QUAD> {
        let mut random_mod_p: Uint<LIMBS_SINGLE>;
        let mut random_mod_q: Uint<LIMBS_SINGLE>;

        loop {
            random_mod_p = Uint::<LIMBS_SINGLE>::random_mod(&mut rng, &NonZero::new(p).unwrap());
            if random_mod_p != q {
                break;
            }
        }

        loop {
            random_mod_q = Uint::<LIMBS_SINGLE>::random_mod(&mut rng, &NonZero::new(q).unwrap());
            if random_mod_q != p {
                break;
            }
        }

        let x = random_mod_p.widening_mul::<LIMBS_SINGLE, LIMBS_DOUBLE>(&random_mod_q);
        x.concat::<LIMBS_QUAD>(&Uint::<LIMBS_DOUBLE>::ZERO)
    }

    /// Generate the full key for the Paillier encryption
    ///
    /// # Parameters
    /// * `rng`: a random number generator for look for primes $p$ and $q$
    ///
    /// # Return
    /// * Secret key ($p$,$q$)
    /// * Public key ($n$,$g$)
    ///
    /// * $p$ large random integer with bit length `LIMBS_SINGLE`
    /// * $q$ large random integer with bit length `LIMBS_SINGLE`
    /// * $n \in \mathbb{Z}_{n^2}^{*}$ public key $p * q$ with bit length `LIMBS_QUAD`
    /// * $g \in \mathcal{B}$ with bit length `LIMBS_QUAD`
    ///
    pub fn generate_key(
        t: u32,
        mut rng: &mut R,
    ) -> KeyPair<LIMBS_SINGLE, LIMBS_DOUBLE, LIMBS_QUAD> {
        let mut p: Uint<LIMBS_SINGLE>;
        let mut q: Uint<LIMBS_SINGLE>;

        loop {
            p = Self::generate_probable_prime(t, &mut rng);
            q = Self::generate_probable_prime(t, &mut rng);

            let (_q1, r1) = p.div_rem(&q.to_nz().unwrap());
            let (_q2, r2) = q.div_rem(&p.to_nz().unwrap());

            if r1.is_zero().unwrap_u8() == 0u8 && r2.is_zero().unwrap_u8() == 0u8 {
                break;
            }
        }

        let n = p.concat::<LIMBS_DOUBLE>(&Uint::ZERO) * q.concat::<LIMBS_DOUBLE>(&Uint::ZERO);

        KeyPair {
            secret_key: SecretKey { p, q },
            public_key: PublicKey { n },
        }
    }
}

#[cfg(test)]
mod tests;
