"""Check research navigation and the actual synchronization scope, not mathematical truth."""
from pathlib import Path
import re
import subprocess
import sys
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]
ALLOWED = {
    '.gitignore', 'README.md', 'README.zh-CN.md', 'RESEARCH_HANDOFF.md', 'AGENTS.md',
    'requirements.txt', 'requirements-audit.txt',
    'docs/01-overview-zh.md', 'docs/02-foundations-zh.md', 'docs/03-results-zh.md',
    'docs/04-evidence-zh.md', 'docs/05-design-zh.md', 'docs/06-process-zh.md',
    'docs/appendices/legacy-models-zh.md', 'docs/appendices/mechanisms-zh.md',
    'tests/check_identities.py', 'tests/check_cut_obstruction.py',
    'tests/check_repository.py', 'tests/check_repository_negative.py',
    'references/README.md', 'audits/sabra/README.md', 'audits/sabra/check_nonlinearity.py',
    'audits/logarithmic/check_review_claims.py',
}
# These are optional local inputs, never synchronization permissions.
LOCAL_MANUSCRIPTS = {
    'manuscript/model.tex', 'manuscript/sparse-modification.tex',
    'manuscript/sparse-modification-zh.tex', 'submission/jsp/README.md',
    'submission/jsp/cover-letter-en.tex', 'submission/jsp/cover-letter-zh.tex',
    'submission/jsp/metadata-en.md', 'submission/jsp/metadata-zh.md',
}
PREVIOUSLY_TRACKED_LOCAL = {'manuscript/model.tex'}
LEGACY = {'docs/'+name for name in (
    '01-model-zh.md', '02-proof-zh.md', '03-novelty-zh.md', '04-development-zh.md',
    '05-validation-zh.md', '06-organization-zh.md', '07-mechanism-screening-zh.md',
    '08-coupled-cascade-zh.md', '09-scale-dependent-exchange-zh.md',
    '10-logarithmic-exchange-zh.md', '11-design-record-zh.md', '12-conversation-log-zh.md',
    '13-proof-dossier-en.md', '14-proof-dossier-zh.md', '15-goal-status-zh.md',
    '16-review-guide-zh.md', '17-review-guide-en.md', '18-independent-review-zh.md',
)}

def require(condition, message):
    if not condition:
        raise AssertionError(message)

def git_names(root, args):
    result = subprocess.run(['git', '-C', str(root), *args, '-z'],
                            capture_output=True, check=True)
    return set(result.stdout.decode('utf-8').split('\0'))-{''}

def anchors(text):
    """GitHub-style headings for the simple headings used in these documents."""
    result, titles = set(), []
    for title in re.findall(r'(?m)^#{1,6} (.+)$', re.sub(r'```.*?```', '', text, flags=re.S)):
        require(title not in titles, 'Duplicate Markdown heading: '+title)
        titles.append(title)
        plain = re.sub(r'\[([^]]+)\]\([^)]*\)', r'\1', title).lower()
        result.add(re.sub(r'[^\w\- ]', '', plain).replace(' ', '-'))
    return result

def read_file(root, name, staged=False):
    path = root/name
    require(path.is_file() and not path.is_symlink(), 'Missing or linked file: '+name)
    if not staged:
        return path.read_text(encoding='utf-8')
    result = subprocess.run(['git', '-C', str(root), 'show', ':'+name],
                            capture_output=True)
    require(result.returncode == 0, 'Missing staged file: '+name)
    return result.stdout.decode('utf-8')

def check_text(name, text):
    if not name.endswith(('.md', '.tex')):
        return
    require(not re.search(r'[\x00-\x08\x0b\x0c\x0e-\x1f]', text), 'Control character: '+name)
    require(not re.search(r'[A-Za-z]:[/\\](?:Users|Program Files|桌面)', text),
            'Machine-local path: '+name)
    if name.endswith('.md'):
        require(len(re.findall(r'(?m)^\$\$\s*$', text)) % 2 == 0,
                'Unpaired display math: '+name)
        require(len(re.findall(r'(?m)^```', text)) % 2 == 0, 'Unpaired code fence: '+name)
        require(not re.search(r'\\m_[wz]\b', text), 'Invalid migrated math macro: '+name)

def check_links(root, contents, heading_ids, permitted):
    for name, text in contents.items():
        if not name.endswith('.md'):
            continue
        prose = re.sub(r'```.*?```', '', text, flags=re.S)
        targets = re.findall(r'\[[^\]]*\]\(([^)]+)\)', prose)
        targets += re.findall(r'(?m)^\[[^\]]+\]:\s*(\S+)', prose)
        for target in targets:
            target = target.strip('<>')
            if re.match(r'https?://|mailto:', target):
                continue
            location, sep, fragment = unquote(target).partition('#')
            resolved = ((root/name).parent/location).resolve() if location else root/name
            require(resolved.is_relative_to(root), 'Link leaves repository: '+name+': '+target)
            rel = resolved.relative_to(root).as_posix()
            require(rel in permitted and resolved.is_file(),
                    'Unpublished or missing link: '+name+': '+target)
            if sep and rel.endswith('.md'):
                require(fragment in heading_ids[rel], 'Unknown heading link: '+name+': '+target)

def check_readme_pair(contents):
    en, zh = contents['README.md'], contents['README.zh-CN.md']
    require('[简体中文](README.zh-CN.md)' in en and '[English](README.md)' in zh,
            'Missing README language switch')
    require('English abstract' not in en and 'English abstract' not in zh,
            'Superseded single README abstract')
    formulas = lambda text: re.findall(r'```math\s*([\s\S]*?)```', text)
    require(formulas(en) == formulas(zh), 'Bilingual README formula mismatch')
    def targets(text):
        return set(re.findall(r'\[[^\]]*\]\(([^)]+)\)', text)) - {'README.md', 'README.zh-CN.md'}
    require(targets(en) == targets(zh), 'Bilingual README navigation mismatch')
    print('Bilingual README formulas, navigation and language switches PASS')

def check_local_manuscripts(root, public_contents, public_headings):
    contents = {name: read_file(root, name) for name in sorted(LOCAL_MANUSCRIPTS)}
    heading_ids = dict(public_headings)
    for name, text in contents.items():
        check_text(name, text)
        if name.endswith('.md'):
            heading_ids[name] = anchors(text)
    check_links(root, contents, heading_ids, ALLOWED | LOCAL_MANUSCRIPTS)
    for name, tex in contents.items():
        if not name.endswith('.tex'):
            continue
        labels = re.findall(r'\\label\{([^}]+)\}', tex)
        require(len(labels) == len(set(labels)), 'Duplicate TeX labels: '+name)
        refs = set(re.findall(r'\\(?:eqref|ref)\{([^}]+)\}', tex))
        require(refs <= set(labels), f'Unknown TeX references: {name}: {refs-set(labels)}')
        bib = set(re.findall(r'\\bibitem(?:\[[^]]*\])?\{([^}]+)\}', tex))
        citations = {item.strip() for group in re.findall(r'\\cite(?:\[[^]]*\])?\{([^}]+)\}', tex) for item in group.split(',')}
        require(citations <= bib, f'Unknown bibliography references: {name}: {citations-bib}')
        environments = []
        for token in re.finditer(r'\\(begin|end)\{([^}]+)\}', tex):
            kind, env = token.groups()
            if kind == 'begin':
                environments.append(env)
            else:
                require(bool(environments) and environments.pop() == env,
                        'Mismatched TeX environment: '+name+': '+env)
        require(not environments, 'Unclosed TeX environment: '+name)
        print(f'{name}: TeX labels: {len(labels)}; bibliography entries: {len(bib)} PASS')

    en = contents['manuscript/sparse-modification.tex']
    zh = contents['manuscript/sparse-modification-zh.tex']
    def display_formulas(tex):
        tex = tex.split(r'\begin{document}', 1)[1].split(r'\end{document}', 1)[0]
        matches = re.findall(r'\\begin\{equation\}([\s\S]*?)\\end\{equation\}|\\\[([\s\S]*?)\\\]', tex)
        return [re.sub(r'\s+', '', re.sub(r'\\text\{[^}]*\}', '', a or b)) for a,b in matches]
    require(display_formulas(en) == display_formulas(zh), 'Bilingual display formula mismatch')
    for command in ('label', 'bibitem', 'cite', 'ref', 'eqref'):
        pattern = r'\\'+command+r'\{([^}]+)\}'
        require(sorted(re.findall(pattern,en)) == sorted(re.findall(pattern,zh)),
                'Bilingual reference mismatch: '+command)
    abstract = re.search(r'\\begin\{abstract\}([\s\S]*?)\\end\{abstract\}', en).group(1)
    abstract_words = len(abstract.split())
    require(150 <= abstract_words <= 250, 'JSP abstract word count: '+str(abstract_words))
    for name in ('manuscript/sparse-modification.tex', 'manuscript/sparse-modification-zh.tex',
                 'submission/jsp/cover-letter-en.tex', 'submission/jsp/cover-letter-zh.tex',
                 'submission/jsp/metadata-en.md', 'submission/jsp/metadata-zh.md'):
        require('Jiang Shirong' in contents[name] and '1296811180@qq.com' in contents[name],
                'Submission author or correspondence missing: '+name)
    for name in ('manuscript/sparse-modification.tex', 'manuscript/sparse-modification-zh.tex'):
        require('not limited to language editing' not in contents[name]
                and '不限于语言润色' not in contents[name], 'Superseded AI wording: '+name)
    print(f'Bilingual display equations: {len(display_formulas(en))}; JSP abstract: {abstract_words} words PASS')

def check(root=ROOT, staged=False, local_manuscripts=False):
    root = root.resolve()
    require(not staged or (root/'.git').exists(), 'Staged check requires a Git repository')
    if (root/'.git').exists():
        indexed = git_names(root, ['ls-files'])
        require(indexed <= ALLOWED, 'Unexpected synchronized files: '+str(sorted(indexed-ALLOWED)))
        if staged:
            entries = subprocess.run(['git', '-C', str(root), 'ls-files', '--stage', '-z'],
                                     capture_output=True, check=True).stdout.decode('utf-8')
            for entry in entries.split('\0'):
                if not entry:
                    continue
                metadata, name = entry.split('\t', 1)
                mode, _, stage = metadata.split()
                require(stage == '0' and mode in {'100644', '100755'},
                        'Linked or unmerged staged file: '+name)
    contents, heading_ids = {}, {}
    for name in sorted(ALLOWED):
        text = read_file(root, name, staged=staged)
        contents[name] = text
        check_text(name, text)
        if name.endswith('.md'):
            heading_ids[name] = anchors(text)
    check_links(root, contents, heading_ids, ALLOWED)
    check_readme_pair(contents)

    if (root/'.git').exists():
        deleted = git_names(root, ['diff', '--cached' if staged else 'HEAD',
                                  '--no-renames', '--name-only', '--diff-filter=D'])
        removable = LEGACY | PREVIOUSLY_TRACKED_LOCAL
        require(deleted <= ALLOWED | removable,
                'Unexpected deletion: '+str(sorted(deleted-ALLOWED-removable)))
        selected = git_names(root, ['diff', '--cached', '--no-renames', '--name-only']) if staged else indexed
        removed = deleted & removable
        require(selected <= ALLOWED | removed,
                'Unexpected synchronized files: '+str(sorted(selected-ALLOWED-removed)))
        if staged:
            require(bool(selected), 'No staged files to check')
        print(f'Git research allowlist: {len(indexed)} PASS; reviewed exclusions: {len(removed)}')
    print(f'Required research files, headings, math fences and Markdown links: {len(ALLOWED)} PASS')
    if local_manuscripts:
        check_local_manuscripts(root, contents, heading_ids)
    else:
        print('Local manuscripts and submission materials are excluded and not required.')
    print('Scope: research files and navigation; not a proof, native TeX compiler, priority or complete secret scanner.')

if __name__ == '__main__':
    check(staged='--staged' in sys.argv, local_manuscripts='--local-manuscripts' in sys.argv)
