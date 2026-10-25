"""Validate the two local bilingual manuscript pairs without publishing their files."""
from check_repository import ROOT, check_local_manuscripts

if __name__ == '__main__':
    check_local_manuscripts(ROOT)
    print('Local source consistency only; native TeX compilation is a separate check.')
