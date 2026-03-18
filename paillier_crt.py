import random
from collections import namedtuple

import gmpy2

"""使用中国剩余定理优化Paillier算法"""
"""Optimize paillier algorithm's use Chinese Remainder Theorem"""


class PaillierKeyGenerator:
    @staticmethod
    def _get_prime_over(N):
        rand_func = random.SystemRandom()
        r = gmpy2.mpz(rand_func.getrandbits(N))
        r = gmpy2.bit_set(r, N - 1)
        return int(gmpy2.next_prime(r))

    @staticmethod
    def _generate_p_q(key_size):
        p = q = None
        n_len = 0
        while n_len != key_size:
            p = PaillierKeyGenerator._get_prime_over(key_size // 2)
            while gmpy2.mod(p, 4) != 3:
                p = PaillierKeyGenerator._get_prime_over(key_size // 2)
            q = p
            while q == p:
                q = PaillierKeyGenerator._get_prime_over(key_size // 2)
                while gmpy2.mod(q, 4) != 3:
                    q = PaillierKeyGenerator._get_prime_over(key_size // 2)
            n = p * q
            n_len = n.bit_length()
        return p, q

    @staticmethod
    def generate_keypair(key_size, s=1):
        p, q = PaillierKeyGenerator._generate_p_q(key_size)

        n = p * q
        lam = (p - 1) * (q - 1) // 2
        mu = gmpy2.invert(lam, n)
        x = random.randint(2, min(p, q))
        h = -pow(x, 2)

        n_square = pow(n, 2)
        hs = gmpy2.powmod(h, n, n_square)
        max_alpha = pow(2, int(gmpy2.ceil(key_size // 2)))

        PublicKey = namedtuple("PublicKey", "n hs max_alpha")
        PrivateKey = namedtuple("PrivateKey", "public_key p q")
        public_key = PublicKey(n=n, hs=hs, max_alpha=max_alpha)
        private_key = PrivateKey(public_key=public_key, p=p, q=q)
        return public_key, private_key


class Paillier:
    CIPHER_MODE_ENCRYPT = 0
    CIPHER_MODE_DECRYPT = 1

    def __init__(self, cipher_mode, cipher_key):
        if cipher_mode == Paillier.CIPHER_MODE_ENCRYPT:
            self.public_key = cipher_key
            self.private_key = None
        elif cipher_mode == Paillier.CIPHER_MODE_DECRYPT:
            self.public_key = cipher_key.public_key
            self.private_key = cipher_key

            self.hp = gmpy2.invert((self.private_key.p - 1) * self.private_key.q, self.private_key.p)
            self.hq = gmpy2.invert((self.private_key.q - 1) * self.private_key.p, self.private_key.q)
            self.p_square = pow(self.private_key.p, 2)
            self.q_square = pow(self.private_key.q, 2)
            self.p_inverse = gmpy2.invert(self.private_key.p, self.private_key.q)
        else:
            raise ValueError('cipher_mode value must be either CIPHER_MODE_ENCRYPT or CIPHER_MODE_DECRYPT')
        self.cipher_mode = cipher_mode
        self.n_square = pow(self.public_key.n, 2)

    def fn_L(self, x, denominator):
        return (x - 1) // denominator

    def encrypt(self, m):
        alpha = random.randint(2, self.public_key.max_alpha - 1)
        # 使用powmod优化模幂运算
        cipher_text = gmpy2.mod((m * self.public_key.n + 1) * gmpy2.powmod(self.public_key.hs, alpha, self.n_square),
                                self.n_square)
        return CryptoNumber(cipher_text, self.n_square)

    def decrypt(self, crypto_number):
        # 使用powmod优化模幂运算
        mp = gmpy2.mod(self.fn_L(gmpy2.powmod(crypto_number.cipher_text, self.private_key.p - 1, self.p_square),
                                 self.private_key.p) * self.hp, self.private_key.p)
        mq = gmpy2.mod(self.fn_L(gmpy2.powmod(crypto_number.cipher_text, self.private_key.q - 1, self.q_square),
                                 self.private_key.q) * self.hq, self.private_key.q)
        return gmpy2.mod(self.crt(mp, mq), self.public_key.n)

    def crt(self, mp, mq):
        u = gmpy2.mod(gmpy2.mul(mq - mp, self.p_inverse), self.private_key.q)
        return mp + (u * self.private_key.p)


class CryptoNumber:
    def __init__(self, cipher_text, n_square):
        self.cipher_text = cipher_text
        self.n_square = n_square

    def __add__(self, other):
        if isinstance(other, CryptoNumber):
            sum_ciphertext = gmpy2.mod(self.cipher_text * other.cipher_text, self.n_square)
            return CryptoNumber(sum_ciphertext, self.n_square)
        else:
            pass

    def __mul__(self, other):
        if isinstance(other, CryptoNumber):
            raise NotImplementedError('not supported between instance of "CryptoNumber" and "CryptoNumber"')
        else:
            mul_cipher_text = gmpy2.mod(pow(self.cipher_text, other), self.n_square)
            return CryptoNumber(mul_cipher_text, self.n_square)


import time

if __name__ == "__main__":
    # 1. 生成密钥对（此处不测量密钥生成时间，单独计时）
    print("正在生成密钥对（2048位）...")
    start_time = time.perf_counter()
    public_key, private_key = PaillierKeyGenerator.generate_keypair(2048)
    keygen_time = time.perf_counter() - start_time
    print(f"密钥生成耗时: {keygen_time:.6f} 秒")
    print(f"模数 n 的实际位数: {public_key.n.bit_length()} bit\n")



    # 2. 创建加密器和解密器
    encryptor = Paillier(Paillier.CIPHER_MODE_ENCRYPT, public_key)
    decryptor = Paillier(Paillier.CIPHER_MODE_DECRYPT, private_key)

    # 3. 测试单个数字的加密解密并计时
    plaintext = 123456789
    print(f"原始明文: {plaintext}")

    # 加密计时
    start_time = time.perf_counter()
    cipher = encryptor.encrypt(plaintext)
    encrypt_time = time.perf_counter() - start_time
    print(f"加密耗时: {encrypt_time:.9f} 秒")
    print(f"加密后密文: {cipher.cipher_text}\n")

    # 解密计时
    start_time = time.perf_counter()
    decrypted = decryptor.decrypt(cipher)
    decrypt_time = time.perf_counter() - start_time
    print(f"解密耗时: {decrypt_time:.9f} 秒")
    print(f"解密结果: {decrypted}\n")

    # 4. 准备两个密文用于同态运算测试
    a, b = 50, 30
    ca = encryptor.encrypt(a)
    cb = encryptor.encrypt(b)

    # 5. 同态加法（密文乘法）计时
    start_time = time.perf_counter()
    c_sum = ca + cb   # 对应密文相乘
    add_time = time.perf_counter() - start_time
    decrypted_sum = decryptor.decrypt(c_sum)
    print(f"{a} + {b} = {decrypted_sum} (同态加法结果)")
    print(f"同态加法（乘法运算）耗时: {add_time:.9f} 秒\n")

    # 6. 标量乘法（密文幂运算）计时
    scalar = 4
    start_time = time.perf_counter()
    c_mul = ca * scalar   # 对应密文的 scalar 次幂
    mul_pow_time = time.perf_counter() - start_time
    decrypted_mul = decryptor.decrypt(c_mul)
    print(f"{a} * {scalar} = {decrypted_mul} (同态乘法结果)")
    print(f"同态乘法（幂运算）耗时: {mul_pow_time:.9f} 秒")

    # 可选：再测试一次模幂运算单独计时（比如单独计算 hs^alpha）
    alpha = random.randint(2, public_key.max_alpha - 1)
    start_time = time.perf_counter()
    _ = gmpy2.powmod(public_key.hs, alpha, encryptor.n_square)
    powmod_time = time.perf_counter() - start_time
    print(f"\n单次模幂运算 (hs^alpha) 耗时: {powmod_time:.9f} 秒")