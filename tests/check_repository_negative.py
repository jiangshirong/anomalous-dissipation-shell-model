"""Exercise upload-scope failures in temporary repositories; leave real files untouched."""
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory
from check_repository import ALLOWED, LOCAL_MANUSCRIPTS, ROOT, check

def stage(root, *names, force=False):
    subprocess.run(['git', '-C', str(root), 'add', *(['-f'] if force else []),
                    '--', *names], check=True, capture_output=True)

def fixture(root, local=False):
    for name in ALLOWED | (LOCAL_MANUSCRIPTS if local else set()):
        target = root/name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT/name, target)
    subprocess.run(['git', '-C', str(root), 'init', '-q'], check=True)
    stage(root, *sorted(ALLOWED))

def append(root, name, text, sync=True):
    path = root/name
    path.write_text(path.read_text(encoding='utf-8')+'\n'+text+'\n', encoding='utf-8')
    if sync:
        stage(root, name)

def verify(root, scenario, expected=None, local=False):
    try:
        with redirect_stdout(StringIO()):
            check(root, staged=True, local_manuscripts=local)
    except AssertionError as error:
        assert expected is not None and expected in str(error), (scenario, str(error))
        print(scenario+': expected failure PASS')
    else:
        assert expected is None, scenario+' was not detected'
        print(scenario+': PASS')

def main():
    cases = {
        'research-only-baseline': None,
        'missing-file': 'Missing or linked file',
        'bad-link': 'Unpublished or missing link',
        'bad-anchor': 'Unknown heading link',
        'unexpected-file': 'Unexpected synchronized files',
        'forced-manuscript': 'Unexpected synchronized files',
        'forced-submission': 'Unexpected synchronized files',
        'forced-pdf': 'Unexpected synchronized files',
        'bad-fence': 'Unpaired code fence',
        'missing-staged-file': 'Missing staged file',
        'staged-bad-link-safe-working-copy': 'Unpublished or missing link',
        'bilingual-readme-formula': 'Bilingual README formula mismatch',
        'bilingual-readme-navigation': 'Bilingual README navigation mismatch',
        'missing-language-switch': 'Missing README language switch',
    }
    for scenario, expected in cases.items():
        with TemporaryDirectory(prefix='shell-research-check-') as folder:
            root = Path(folder)
            fixture(root)
            if scenario == 'missing-file':
                (root/'docs/03-results-zh.md').unlink()
            elif scenario in {'bad-link', 'bad-anchor', 'staged-bad-link-safe-working-copy'}:
                path = root/'README.md'
                original = path.read_text(encoding='utf-8')
                target = 'docs/03-results-zh.md#missing-heading' if scenario == 'bad-anchor' else 'docs/missing.md'
                append(root, 'README.md', '[negative test]('+target+')')
                if scenario == 'staged-bad-link-safe-working-copy':
                    path.write_text(original, encoding='utf-8')
            elif scenario in {'unexpected-file', 'forced-manuscript', 'forced-submission', 'forced-pdf'}:
                name = {'unexpected-file': 'unapproved.txt',
                        'forced-manuscript': 'manuscript/model.tex',
                        'forced-submission': 'submission/jsp/metadata-en.md',
                        'forced-pdf': 'output/pdf/paper.pdf'}[scenario]
                path = root/name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text('temporary fixture only\n', encoding='utf-8')
                stage(root, name, force=True)
            elif scenario == 'bad-fence':
                append(root, 'README.md', '```')
            elif scenario == 'missing-staged-file':
                subprocess.run(['git', '-C', str(root), 'rm', '--cached', '--', 'README.zh-CN.md'],
                               check=True, capture_output=True)
            elif scenario == 'bilingual-readme-formula':
                append(root, 'README.zh-CN.md', '```math\n1=2\n```')
            elif scenario == 'bilingual-readme-navigation':
                append(root, 'README.md', '[dependencies](requirements.txt)')
            elif scenario == 'missing-language-switch':
                path = root/'README.zh-CN.md'
                path.write_text(path.read_text(encoding='utf-8').replace('[English](README.md)', 'English'),
                                encoding='utf-8')
                stage(root, 'README.zh-CN.md')
            verify(root, scenario, expected)

    local_cases = {
        'local-manuscripts-baseline': None,
        'new-tex-ref': 'Unknown TeX references',
        'new-tex-cite': 'Unknown bibliography references',
        'missing-translation': 'Missing or linked file',
        'bilingual-formula': 'Bilingual display formula mismatch',
        'bad-tex-environment': 'Unclosed TeX environment',
    }
    local_available = all((ROOT/name).is_file() for name in LOCAL_MANUSCRIPTS)
    if '--local-manuscripts' in sys.argv:
        assert local_available, 'Local manuscript negative checks require the excluded local materials'
    if not local_available:
        local_cases = {}
        print('Local manuscript negative checks skipped: excluded materials are absent.')
    for scenario, expected in local_cases.items():
        with TemporaryDirectory(prefix='shell-local-check-') as folder:
            root = Path(folder)
            fixture(root, local=True)
            if scenario in {'new-tex-ref', 'new-tex-cite'}:
                invalid = r'\eqref{missing-label}' if scenario == 'new-tex-ref' else r'\cite{missing-source}'
                append(root, 'manuscript/sparse-modification.tex', invalid, sync=False)
            elif scenario == 'missing-translation':
                (root/'submission/jsp/cover-letter-zh.tex').unlink()
            elif scenario == 'bilingual-formula':
                path = root/'manuscript/sparse-modification-zh.tex'
                path.write_text(path.read_text(encoding='utf-8').replace(
                    r'\end{document}', r'\[1=2\]'+'\n'+r'\end{document}', 1), encoding='utf-8')
            elif scenario == 'bad-tex-environment':
                append(root, 'submission/jsp/cover-letter-en.tex', r'\begin{proof}', sync=False)
            verify(root, scenario, expected, local=True)
    baselines = 1 + bool(local_cases)
    print(f'Negative scenarios: {len(cases)+len(local_cases)-baselines} PASS; {baselines} baselines PASS; real workspace untouched.')

if __name__ == '__main__':
    main()
