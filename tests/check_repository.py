"""Validate reviewed publication files and navigation, not mathematical truth."""
from pathlib import Path, PurePosixPath
import argparse
import json
import posixpath
import re
import subprocess
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = 'tests/public-files.json'
FORBIDDEN_PARTS = {'local', '.local-backup', '.git', 'tmp', 'manuscript', 'submission', 'output', '__pycache__'}
FORBIDDEN_SUFFIXES = {'.pdf', '.tex', '.bib', '.bbl', '.cls', '.sty', '.zip', '.bundle', '.aux', '.log', '.gz', '.tar', '.pyc'}
PRIVATE_FILENAMES = {'.env', 'id_rsa', 'id_ed25519', '.npmrc', '.pypirc', 'credentials.json', 'secrets.json'}
ESSENTIAL = {MANIFEST, '.gitignore', 'README.md', 'README.zh-CN.md', 'AGENTS.md', 'RESEARCH_HANDOFF.md'}


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def git(root, *args):
    return subprocess.run(['git', '-C', str(root), *args], capture_output=True, check=True).stdout


def read_file(root, name, staged=False):
    """Staged validation reads index bytes, even when the work tree looks correct."""
    if staged:
        result = subprocess.run(['git', '-C', str(root), 'show', ':' + name], capture_output=True)
        require(result.returncode == 0, 'Missing staged file: ' + name)
        data = result.stdout
    else:
        p = root / name
        require(p.is_file() and not p.is_symlink(), 'Missing or linked file: ' + name)
        require(p.resolve().is_relative_to(root.resolve()), 'File leaves repository: ' + name)
        data = p.read_bytes()
    try:
        return data.decode('utf-8-sig')
    except UnicodeDecodeError as error:
        raise AssertionError('Public file is not UTF-8 text: ' + name) from error


def manifest_files(root=ROOT, staged=False):
    try:
        record = json.loads(read_file(root, MANIFEST, staged))
    except json.JSONDecodeError as error:
        raise AssertionError('Invalid publication manifest JSON') from error
    require(isinstance(record, dict) and set(record) == {'files'}, 'Invalid publication manifest shape')
    names = record['files']
    require(isinstance(names, list) and all(isinstance(n, str) for n in names), 'Invalid publication file list')
    require(len(names) == len(set(names)), 'Duplicate publication manifest path')
    for name in names:
        p = PurePosixPath(name)
        require(name and '\\' not in name and not p.is_absolute() and p.as_posix() == name
                and ':' not in name and all(x not in {'', '.', '..'} for x in p.parts),
                'Noncanonical publication path: ' + name)
        require(not any(part.lower() in FORBIDDEN_PARTS for part in p.parts)
                and p.suffix.lower() not in FORBIDDEN_SUFFIXES
                and p.name.lower() not in PRIVATE_FILENAMES and not p.name.lower().startswith('.env.'),
                'Forbidden publication path: ' + name)
    allowed = set(names)
    require(ESSENTIAL <= allowed, 'Missing essential publication files: ' + str(sorted(ESSENTIAL - allowed)))
    return allowed


def without_fences(text):
    output = []
    fence = None
    for line in text.splitlines():
        match = re.match(r'^ {0,3}(`{3,}|~{3,})(.*)$', line)
        if match:
            token = match.group(1)
            if fence is None:
                fence = token
            elif token[0] == fence[0] and len(token) >= len(fence) and not match.group(2).strip():
                fence = None
            continue
        if fence is None:
            output.append(line)
    require(fence is None, 'Unpaired code fence')
    return '\n'.join(output)


def anchors(text):
    """Heading fragments, including GitHub duplicate-heading numeric suffixes."""
    clean = without_fences(text)
    ids = set(re.findall(r'<a\s+[^>]*(?:id|name)=["\']([^"\']+)["\']', clean, re.I))
    used = set(ids)
    for title in re.findall(r'(?m)^ {0,3}#{1,6}\s+(.+?)\s*#*\s*$', clean):
        title = re.sub(r'\[([^]]+)\]\([^)]*\)', r'\1', title)
        title = re.sub(r'<[^>]*>', '', title)
        title = title.replace('`', '').replace('*', '')
        slug = re.sub(r'[^\w\- ]', '', title.lower()).replace(' ', '-')
        candidate = slug
        suffix = 1
        while candidate in used:
            candidate = slug + '-' + str(suffix)
            suffix += 1
        used.add(candidate)
        ids.add(candidate)
    return ids


def prose_without_math(text):
    """Keep ordinary/code prose while removing declared Markdown/TeX math."""
    lines = []
    math_fence = None
    for line in text.splitlines():
        if math_fence is None:
            match = re.match(r'^ {0,3}(`{3,}|~{3,})math\s*$', line)
            if match:
                math_fence = match.group(1)
                continue
            lines.append(line)
        elif re.match(r'^ {0,3}' + re.escape(math_fence[0]) + '{' + str(len(math_fence)) + r',}\s*$', line):
            math_fence = None
    result = '\n'.join(lines)
    result = re.sub(r'(?<!\\)\$\$[\s\S]*?(?<!\\)\$\$', '', result)
    result = re.sub(r'(?<!\\)\$(?!\$)[^\n]*?(?<!\\)\$', '', result)
    result = re.sub(r'\\\[[\s\S]*?\\\]', '', result)
    result = re.sub(r'\\\([^\n]*?\\\)', '', result)
    return result


def check_text(name, text):
    require(not re.search(r'[\x00-\x08\x0b\x0c\x0e-\x1f]', text), 'Control character: ' + name)
    if name.endswith(('.md', '.json', '.txt')):
        known_root = r'(?i)(?<![A-Za-z0-9])[A-Z]:[/\\](?:Users|Program Files|Windows|桌面|桌面文件夹)(?:[/\\]|\b)'
        drive_path = r'(?i)(?<![A-Za-z0-9])[A-Z]:[/\\][A-Za-z0-9_\u0080-\uffff]'
        require(not re.search(known_root, text) and not re.search(drive_path, prose_without_math(text)),
                'Machine-local path: ' + name)
    credential_patterns = [r'gh[pousr]_[A-Za-z0-9]{20,}', r'github_pat_[A-Za-z0-9_]{20,}',
                           r'sk-[A-Za-z0-9]{20,}', r'-----BEGIN (?:[A-Z ]+ )?PRIVATE KEY-----']
    require(not any(re.search(pattern, text) for pattern in credential_patterns), 'Credential-like material: ' + name)
    if name.endswith('.md'):
        prose = without_fences(text)
        require(len(re.findall(r'(?m)^\s*\$\$\s*$', prose)) % 2 == 0, 'Unpaired display math: ' + name)


def links(text):
    prose = without_fences(text)
    inline = re.findall(r'!?\[[^\]]*\]\((<[^>]*>|[^\s)]+)(?:\s+["\'][^"\']*["\'])?\)', prose)
    reference = re.findall(r'(?m)^ {0,3}\[(?!\^)[^\]]+\]:\s*(\S+)', prose)
    return {target.strip('<>') for target in inline + reference}


def check_links(contents, heading_ids, permitted):
    directories = {parent.as_posix() for name in permitted for parent in PurePosixPath(name).parents}
    for name, text in contents.items():
        if not name.endswith('.md'):
            continue
        for target in links(text):
            decoded = unquote(target)
            if re.match(r'^(?:https?://|mailto:)', decoded, re.I):
                continue
            require(not urlsplit(decoded).scheme and not decoded.startswith('//'), 'Unsupported link scheme: ' + name + ': ' + target)
            location, separator, fragment = decoded.partition('#')
            require(not location.startswith('/') and '\\' not in location, 'Link leaves repository: ' + name + ': ' + target)
            resolved = posixpath.normpath(posixpath.join(posixpath.dirname(name), location)) if location else name
            require(resolved != '..' and not resolved.startswith('../'), 'Link leaves repository: ' + name + ': ' + target)
            require(resolved in permitted or resolved in directories, 'Unpublished or missing link: ' + name + ': ' + target)
            if separator and fragment:
                require(resolved in heading_ids and fragment in heading_ids[resolved], 'Unknown heading link: ' + name + ': ' + target)


def check_readme_pair(contents):
    en, zh = contents['README.md'], contents['README.zh-CN.md']
    require('[简体中文](README.zh-CN.md)' in en and '[English](README.md)' in zh, 'Missing README language switch')
    require('English abstract' not in en and 'English abstract' not in zh, 'Superseded single README abstract')
    formulas = lambda text: re.findall(r'```math\s*([\s\S]*?)```', text)
    require(formulas(en) == formulas(zh), 'Bilingual README formula mismatch')
    navigation = lambda text: links(text) - {'README.md', 'README.zh-CN.md'}
    require(navigation(en) == navigation(zh), 'Bilingual README navigation mismatch')


def check_ignore(text, allowed):
    rules = [line.strip() for line in text.splitlines() if line.strip() and not line.lstrip().startswith('#')]
    require(rules and rules[0] == '*', 'Ignore policy must start with deny-all rule')
    ancestors = {parent.as_posix() for name in allowed for parent in PurePosixPath(name).parents if parent.as_posix() != '.'}
    expected = {'*'} | {'!' + name for name in allowed} | {'!' + parent + '/' for parent in ancestors}
    require(set(rules) == expected and len(rules) == len(set(rules)), 'Ignore policy does not match explicit public allowlist')


def check_index(root, allowed):
    records = git(root, 'ls-files', '--stage', '-z').decode('utf-8').split('\0')
    indexed = set()
    for record in filter(None, records):
        metadata, name = record.split('\t', 1)
        mode, _blob, stage = metadata.split()
        require(stage == '0', 'Unmerged index file: ' + name)
        require(mode in {'100644', '100755'}, 'Linked or unsupported index file: ' + name)
        indexed.add(name)
    require(indexed == allowed, 'Index differs from publication manifest: extra=' + str(sorted(indexed - allowed))
            + '; missing=' + str(sorted(allowed - indexed)))


def check_local_manuscripts(root):
    pairs = [
        ('sparse-modification', 'sparse-modification.tex', 'sparse-modification-zh.tex'),
        ('correlated-environments', 'd1-en.tex', 'd1-zh.tex'),
    ]
    for folder, en_name, zh_name in pairs:
        base = root / 'local/papers' / folder / 'sources'
        texts = []
        for name in (en_name, zh_name):
            path = base / name
            require(path.is_file() and not path.is_symlink(), 'Missing local manuscript: ' + str(path.relative_to(root)))
            text = path.read_text(encoding='utf-8')
            labels = re.findall(r'\\label\{([^}]+)\}', text)
            require(len(labels) == len(set(labels)), 'Duplicate local TeX label: ' + name)
            refs = set(re.findall(r'\\(?:eqref|ref)\{([^}]+)\}', text))
            require(refs <= set(labels), 'Unknown local TeX reference: ' + name)
            texts.append(text)
        equations = lambda text: [re.sub(r'\s+', '', re.sub(r'\\text\{[^}]*\}', '', value)) for value in re.findall(r'\\begin\{equation\}([\s\S]*?)\\end\{equation\}', text)]
        require(equations(texts[0]) == equations(texts[1]), 'Bilingual local numbered-equation mismatch: ' + folder)
        print('Local bilingual manuscript pair PASS: ' + folder)


def check(root=ROOT, staged=False, local_manuscripts=False):
    root = Path(root).resolve()
    allowed = manifest_files(root, staged)
    if staged:
        check_index(root, allowed)
    contents = {name: read_file(root, name, staged) for name in sorted(allowed)}
    for name, text in contents.items():
        check_text(name, text)
    heading_ids = {name: anchors(text) for name, text in contents.items() if name.endswith('.md')}
    check_links(contents, heading_ids, allowed)
    check_readme_pair(contents)
    check_ignore(contents['.gitignore'], allowed)
    if not staged and (root / '.git').exists():
        for name in allowed:
            result = subprocess.run(['git', '-C', str(root), 'check-ignore', '--no-index', '-q', '--', name], capture_output=True)
            require(result.returncode == 1, 'Approved public file is ignored: ' + name)
        for name in ('local/papers/example.pdf', 'manuscript/example.tex', 'submission/example.txt', 'tmp/example.py', 'unreviewed.txt', 'docs/unreviewed.md'):
            result = subprocess.run(['git', '-C', str(root), 'check-ignore', '--no-index', '-q', '--', name], capture_output=True)
            require(result.returncode == 0, 'Private or unreviewed file is not ignored: ' + name)
    if local_manuscripts:
        check_local_manuscripts(root)
    print('Publication manifest, navigation, bilingual README and ignore policy PASS: ' + str(len(allowed)) + ' files')
    print('Mode: ' + ('actual Git index' if staged else 'working tree; historical staged deletions are not required'))
    print('Scope: publication structure and finite consistency; not an infinite-dimensional proof, peer review, priority certification or exhaustive secret scan.')
    return allowed


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--staged', action='store_true')
    parser.add_argument('--local-manuscripts', action='store_true')
    args = parser.parse_args()
    check(staged=args.staged, local_manuscripts=args.local_manuscripts)
