# -*- coding: utf-8 -*-
"""각 주차의 마크다운을 그 주차 폴더의 guides.js 에 내장한다.
file:// 로 열면 브라우저가 fetch 를 막으므로, 안내문을 스크립트로 들고 있어야 한다.
마크다운을 고친 뒤에는 이 스크립트를 다시 돌린다.

    python tools/build-guides.py
"""
import io, json, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WEEKS = [
    ('week10', ['README.md', 'worksheet.md', 'ncbi-guide.md']),
    ('week11', ['README.md', 'worksheet.md']),
]

for week, files in WEEKS:
    base = os.path.join(ROOT, 'docs', week)
    out = []
    for rel in files:
        p = os.path.join(base, *rel.split('/'))
        if not os.path.exists(p):
            sys.exit('없는 파일: ' + p)
        out.append((rel, io.open(p, encoding='utf-8').read()))

    body = ['// 자동 생성 파일 — 직접 고치지 말 것. tools/build-guides.py 로 다시 만든다.',
            '// 이 주차의 마크다운을 내장해 file:// 에서도 안내문이 열리게 한다.',
            'window.GUIDES = {']
    body += ['  %s: %s,' % (json.dumps(k, ensure_ascii=False), json.dumps(v, ensure_ascii=False))
             for k, v in out]
    body += ['};', '']

    dst = os.path.join(base, 'guides.js')
    io.open(dst, 'w', encoding='utf-8', newline='\n').write('\n'.join(body))
    print('%s/guides.js  %d bytes · %s' % (week, os.path.getsize(dst),
          ', '.join('%s %d자' % (k, len(v)) for k, v in out)))
