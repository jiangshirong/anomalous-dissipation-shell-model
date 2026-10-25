"""Exercise publication-scope failures in disposable repositories only."""
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
import argparse
import json
import shutil
import subprocess
from check_repository import ROOT, MANIFEST, check, manifest_files, anchors, links


def command(root, *args, input_data=None):
    return subprocess.run(['git', '-C', str(root), *args], input=input_data, capture_output=True, check=True).stdout


def stage(root, *names, force=False):
    command(root, 'add', *(['-f'] if force else []), '--', *names)


def fixture(root):
    allowed = manifest_files(ROOT)
    for name in allowed:
        target = root / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / name, target)
    command(root, 'init', '-q')
    command(root, 'config', 'core.autocrlf', 'false')
    stage(root, *sorted(allowed))
    return allowed


def append(root, name, content, staged=True):
    path = root / name
    path.write_text(path.read_text(encoding='utf-8-sig') + '\n' + content + '\n', encoding='utf-8')
    if staged:
        stage(root, name)


def verify(root, scenario, expected=None, staged=True):
    try:
        with redirect_stdout(StringIO()):
            check(root, staged=staged)
    except AssertionError as error:
        assert expected is not None and expected in str(error), (scenario, str(error), expected)
        print(scenario + ': expected failure PASS')
    else:
        assert expected is None, scenario + ' was not detected'
        print(scenario + ': PASS')


def main():
    cases = {
        'staged-baseline': None,
        'working-tree-baseline': None,
        'legitimate-mathematical-colons': None,
        'missing-file': 'Missing or linked file',
        'bad-link': 'Unpublished or missing link',
        'bad-anchor': 'Unknown heading link',
        'outside-repository-link': 'Link leaves repository',
        'forced-local-paper': 'Index differs from publication manifest',
        'forced-unapproved-file': 'Index differs from publication manifest',
        'private-path-in-manifest': 'Forbidden publication path',
        'paper-path-in-manifest': 'Forbidden publication path',
        'bad-fence': 'Unpaired code fence',
        'bad-display-math': 'Unpaired display math',
        'missing-language-switch': 'Missing README language switch',
        'bilingual-formula': 'Bilingual README formula mismatch',
        'bilingual-navigation': 'Bilingual README navigation mismatch',
        'staged-bad-working-good': 'Unpublished or missing link',
        'staged-ignore-bad-working-good': 'Ignore policy does not match',
        'missing-staged-file': 'Index differs from publication manifest',
        'unmerged-file': 'Unmerged index file',
        'index-symlink': 'Linked or unsupported index file',
        'machine-local-path': 'Machine-local path',
        'machine-local-backslash-path': 'Machine-local path',
        'machine-local-path-in-math': 'Machine-local path',
        'credential-pattern': 'Credential-like material',
        'broad-ignore-exception': 'Ignore policy does not match',
        'noncanonical-manifest-path': 'Noncanonical publication path',
    }
    for scenario, expected in cases.items():
        with TemporaryDirectory(prefix='shell-public-check-') as folder:
            root = Path(folder)
            fixture(root)
            working = scenario in {'working-tree-baseline', 'missing-file'}
            if scenario == 'missing-file':
                (root / 'README.md').unlink()
            elif scenario in {'bad-link', 'bad-anchor', 'outside-repository-link', 'staged-bad-working-good'}:
                original = (root / 'README.md').read_bytes()
                target = {'bad-anchor': 'README.md#missing-section', 'outside-repository-link': '../private.md'}.get(scenario, 'docs/missing-file.md')
                append(root, 'README.md', '[negative fixture](' + target + ')')
                if scenario == 'staged-bad-working-good':
                    (root / 'README.md').write_bytes(original)
            elif scenario in {'forced-local-paper', 'forced-unapproved-file'}:
                name = 'local/papers/private.pdf' if scenario == 'forced-local-paper' else 'unreviewed.txt'
                path = root / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text('Disposable fixture only\n', encoding='utf-8')
                stage(root, name, force=True)
            elif scenario in {'private-path-in-manifest', 'paper-path-in-manifest', 'noncanonical-manifest-path'}:
                name = {'private-path-in-manifest': 'local/notes.md', 'paper-path-in-manifest': 'docs/paper.pdf',
                        'noncanonical-manifest-path': 'docs/../private.md'}[scenario]
                path = root / MANIFEST
                obj = json.loads(path.read_text(encoding='utf-8-sig'))
                obj['files'].append(name)
                path.write_text(json.dumps(obj), encoding='utf-8')
                stage(root, MANIFEST)
            elif scenario == 'bad-fence':
                append(root, 'README.md', '```')
            elif scenario == 'bad-display-math':
                append(root, 'README.md', '$$')
            elif scenario == 'missing-language-switch':
                p = root / 'README.zh-CN.md'
                p.write_text(p.read_text(encoding='utf-8-sig').replace('[English](README.md)', 'English'), encoding='utf-8')
                stage(root, 'README.zh-CN.md')
            elif scenario == 'bilingual-formula':
                append(root, 'README.zh-CN.md', '```math\n1=2\n```')
            elif scenario == 'bilingual-navigation':
                title = sorted(anchors((root / 'README.md').read_text(encoding='utf-8-sig')))[0]
                append(root, 'README.md', '[readme title](README.md#' + title + ')')
            elif scenario in {'staged-ignore-bad-working-good', 'broad-ignore-exception'}:
                original = (root / '.gitignore').read_bytes()
                append(root, '.gitignore', '!docs/**')
                if scenario == 'staged-ignore-bad-working-good':
                    (root / '.gitignore').write_bytes(original)
            elif scenario == 'missing-staged-file':
                command(root, 'rm', '--cached', '--', 'README.md')
            elif scenario in {'unmerged-file', 'index-symlink'}:
                if scenario == 'index-symlink':
                    blob = command(root, 'hash-object', '-w', '--stdin', input_data=b'../private\n').decode().strip()
                    command(root, 'update-index', '--cacheinfo', '120000,' + blob + ',README.md')
                else:
                    blob = command(root, 'rev-parse', ':README.md').decode().strip()
                    data = ('0 ' + '0' * 40 + '\tREADME.md\n100644 ' + blob + ' 1\tREADME.md\n').encode()
                    command(root, 'update-index', '--index-info', input_data=data)
            elif scenario == 'legitimate-mathematical-colons':
                append(root, 'README.md', '```math\n' + r'V=\{a:\sum_n4^na_n^2<\infty\},\quad\sum_{j:\,0<\theta_j<1}r_j^{-1}' + '\n```\n' + r'Inline condition `j:\,0<\theta_j<1`.')
                append(root, 'README.zh-CN.md', '```math\n' + r'V=\{a:\sum_n4^na_n^2<\infty\},\quad\sum_{j:\,0<\theta_j<1}r_j^{-1}' + '\n```\n' + r'行内条件 `j:\,0<\theta_j<1`。')
            elif scenario in {'machine-local-path', 'machine-local-backslash-path', 'machine-local-path-in-math'}:
                value = 'C:' + ('/Users/example/private' if scenario != 'machine-local-backslash-path' else '\\Users\\example\\private')
                append(root, 'README.md', '```math\n' + value + '\n```' if scenario == 'machine-local-path-in-math' else value)
            elif scenario == 'credential-pattern':
                append(root, 'README.md', 'gh' + 'p_' + 'A' * 30)
            verify(root, scenario, expected, staged=not working)
    print('PASS: ' + str(sum(value is not None for value in cases.values())) + ' negative scenarios and 3 baselines; real index untouched.')


if __name__ == '__main__':
    main()
