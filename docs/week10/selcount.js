/* 서열을 끌어 고르면 고른 길이를 띄운다.
   SELCOUNT.watch(상자) 로 상자를 등록한다.
   세는 규칙
     - class="num" "ln" "lab" "cut" 안의 글자는 세지 않는다
     - class="aa"  안의 글자는 아미노산으로 따로 센다
     - 그 밖의 A·C·G·T 는 핵산으로 센다 */
(function (root) {
  var boxes = [], badge = null, pending = false;
  // 서열이 아닌 글자: 위치 번호(num·ln), 효소 이름표(lab), 절단 표시(cut)
  var SKIP = ['num', 'ln', 'lab', 'cut'];

  function tagOf(node) {
    for (var e = node.nodeType === 1 ? node : node.parentNode; e && e !== document.body; e = e.parentNode) {
      if (!e.classList) continue;
      for (var i = 0; i < SKIP.length; i++) if (e.classList.contains(SKIP[i])) return 'skip';
      if (e.classList.contains('aa')) return 'aa';
    }
    return 'nt';
  }
  function inBox(node) {
    for (var e = node.nodeType === 1 ? node : node.parentNode; e; e = e.parentNode) {
      if (boxes.indexOf(e) >= 0) return true;
    }
    return false;
  }
  function count(range) {
    var n = { nt: 0, aa: 0 }, host = range.commonAncestorContainer;
    if (host.nodeType === 3) {
      add(n, tagOf(host), range.toString());
      return n;
    }
    var w = document.createTreeWalker(host, NodeFilter.SHOW_TEXT, null), node;
    while ((node = w.nextNode())) {
      if (!range.intersectsNode(node)) continue;
      var s = node.nodeValue;
      var from = node === range.startContainer ? range.startOffset : 0;
      var to = node === range.endContainer ? range.endOffset : s.length;
      if (to > from) add(n, tagOf(node), s.slice(from, to));
    }
    return n;
  }
  function add(n, tag, s) {
    if (tag === 'skip') return;
    if (tag === 'aa') { var a = s.match(/[A-Za-z*]/g); n.aa += a ? a.length : 0; }
    else { var t = s.match(/[ACGTacgt]/g); n.nt += t ? t.length : 0; }
  }
  function ensure() {
    if (badge) return badge;
    badge = document.createElement('div');
    badge.className = 'selcount';
    badge.setAttribute('aria-hidden', 'true');
    document.body.appendChild(badge);
    return badge;
  }
  function hide() { if (badge) badge.style.display = 'none'; }
  function show(range, text) {
    var rects = range.getClientRects();
    var r = rects.length ? rects[rects.length - 1] : range.getBoundingClientRect();
    if (!r || (!r.width && !r.height)) return hide();
    var b = ensure();
    b.textContent = text;
    b.style.display = 'block';
    var w = b.offsetWidth, h = b.offsetHeight;
    var x = Math.min(Math.max(r.right, w / 2 + 6), window.innerWidth - w / 2 - 6);
    var y = r.top - h - 8;
    if (y < 4) y = r.bottom + 8;          // 위가 좁으면 아래에 붙인다
    b.style.left = Math.round(x) + 'px';
    b.style.top = Math.round(y) + 'px';
  }
  function update() {
    pending = false;
    var sel = window.getSelection();
    if (!sel || !sel.rangeCount || sel.isCollapsed) return hide();
    var range = sel.getRangeAt(0);
    if (!inBox(range.startContainer) || !inBox(range.endContainer)) return hide();
    var n = count(range), parts = [];
    if (n.nt) parts.push(n.nt.toLocaleString() + ' nt');
    if (n.aa) parts.push(n.aa.toLocaleString() + ' aa');
    if (!parts.length) return hide();
    show(range, parts.join(' · '));
  }
  function queue() {
    if (pending) return;
    pending = true;
    // requestAnimationFrame 은 창이 뒤에 있으면 돌지 않는다. 타이머로 미룬다.
    setTimeout(update, 0);
  }

  document.addEventListener('selectionchange', queue);
  window.addEventListener('scroll', queue, true);
  window.addEventListener('resize', queue);

  root.SELCOUNT = {
    watch: function (el) { if (el && boxes.indexOf(el) < 0) boxes.push(el); }
  };
})(window);
