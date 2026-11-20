/* App shell: one page ("lab") per topic. Problems are presets of the lab, chosen from a
   dropdown on the page. Routes: #/  ·  #/lab/<topic>[/<problemId>]  ·  #/sources  ·  #/selftest
   (old #/p/<problemId> links redirect). */
(function () {
  'use strict';
  const { el, tex } = QM.ui;
  const R = QM.registry;
  const main = document.getElementById('main');
  const nav = document.getElementById('nav');
  let mounted = null;

  // ---------- theme ----------
  const themeKey = 'qm-sim-theme';
  try { const t = localStorage.getItem(themeKey); if (t) document.documentElement.dataset.theme = t; } catch (e) { /* ignore */ }
  document.getElementById('themeBtn').addEventListener('click', () => {
    const dark = document.documentElement.dataset.theme ? document.documentElement.dataset.theme === 'dark' : matchMedia('(prefers-color-scheme: dark)').matches;
    document.documentElement.dataset.theme = dark ? 'light' : 'dark';
    try { localStorage.setItem(themeKey, document.documentElement.dataset.theme); } catch (e) { /* ignore */ }
    if (mounted && mounted.redraw) mounted.redraw(); else route();
  });

  const topicProblems = (tid) => R.problems.filter((p) => p.topic === tid);
  const labHref = (p) => `#/lab/${p.topic}/${p.id}`;

  // ---------- sidebar ----------
  function buildNav(filter = '') {
    nav.innerHTML = '';
    const f = filter.trim().toLowerCase();
    nav.append(el('a', { href: '#/' }, 'Overview'));
    nav.append(el('div', { class: 'nav-topic' }, 'Labs'));
    R.topics.forEach((t) => {
      const ps = topicProblems(t.id).filter((p) => !p.sandbox && (!f || (t.title + ' ' + p.title + ' ' + (p.source || '') + ' ' + (p.sheet || '')).toLowerCase().includes(f)));
      if (f && !ps.length && !t.title.toLowerCase().includes(f)) return;
      nav.append(el('a', { href: '#/lab/' + t.id, class: 'nav-lab', 'data-topic': t.id }, t.title));
      ps.forEach((p) => nav.append(el('a', { href: labHref(p), class: 'nav-sub' }, p.title, p.sheet ? el('span', { class: 'src' }, p.sheet) : null)));
    });
    nav.append(el('div', { class: 'nav-topic' }, 'Meta'));
    nav.append(el('a', { href: '#/sources' }, 'Sources & progress'));
    nav.append(el('a', { href: '#/selftest' }, 'Self-test'));
    markActive();
  }
  function markActive() {
    const h = location.hash || '#/';
    const topic = (h.match(/^#\/lab\/([^/]+)/) || [])[1];
    nav.querySelectorAll('a').forEach((a) => {
      const href = a.getAttribute('href');
      a.classList.toggle('active', href === h || (a.classList.contains('nav-lab') && a.dataset.topic === topic && h === '#/lab/' + topic));
      a.classList.toggle('in-lab', a.classList.contains('nav-lab') && a.dataset.topic === topic);
    });
  }
  document.getElementById('search').addEventListener('input', (e) => buildNav(e.target.value));

  function unmount() { if (mounted && mounted.destroy) try { mounted.destroy(); } catch (e) { console.error(e); } mounted = null; }

  // ---------- pages ----------
  function home() {
    main.append(el('div', { class: 'prob-head' },
      el('h2', null, 'QM I – problem simulations'),
      el('div', { class: 'statement', html: 'Each topic is one general <b>lab</b>: choose the Hamiltonian, its parameters, the number of dots/spins and the initial state freely. The problems from <code>QM_I/Problems</code> are <i>presets</i> of these labs – pick one from the dropdown at the top of the lab to load its values, read the statement and reveal the solution. Units: $\\hbar = 1$.' })));
    const cards = el('div', { class: 'cards' });
    R.topics.forEach((t) => {
      const ps = topicProblems(t.id).filter((p) => !p.sandbox);
      cards.append(el('div', { class: 'card' }, el('h3', null, el('a', { href: '#/lab/' + t.id }, t.title)), el('div', { class: 'note', html: t.desc }),
        el('ul', null, ps.map((p) => el('li', null, el('a', { href: labHref(p) }, p.title))))));
    });
    main.append(cards);
  }

  function sourcesPage() {
    main.append(el('div', { class: 'prob-head' }, el('h2', null, 'Sources & progress'),
      el('div', { class: 'statement', html: 'Files in <code>QM_I/Problems</code> and the lab presets built from them. When new files are added, they are compared against <code>simulace/PROCESSED.md</code> and this list.' })));
    const tbl = el('table', { class: 'sources' }, el('tr', null, el('th', null, 'File'), el('th', null, 'Content'), el('th', null, 'Status'), el('th', null, 'Presets')));
    R.sources.forEach((s) => {
      tbl.append(el('tr', null,
        el('td', null, el('a', { href: '../Problems/' + encodeURIComponent(s.file), target: '_blank' }, s.file)),
        el('td', { html: s.title + (s.note ? `<div class="note">${s.note}</div>` : '') }),
        el('td', null, el('span', { class: 'status ' + s.status }, s.status)),
        el('td', null, s.problems.map((id) => { const p = R.problems.find((q) => q.id === id); return el('div', null, p ? el('a', { href: labHref(p) }, p.title) : el('span', { class: 'status todo' }, id + ' (missing)')); }))));
    });
    main.append(tbl);
  }

  function selftest() {
    main.append(el('div', { class: 'prob-head' }, el('h2', null, 'Self-test'), el('div', { class: 'statement' }, 'Numerical checks of the labs against the analytic results of the problems.')));
    const list = el('ol', { class: 'selftest panel', style: { padding: '12px 32px' } });
    main.append(list);
    let pass = 0;
    QM.tests.forEach((t) => {
      let ok = false, detail = '';
      try { const r = t.fn(); ok = r === true || (r && r.pass); detail = r && r.detail ? r.detail : ''; } catch (e) { detail = 'Exception: ' + e.message; console.error(e); }
      if (ok) pass++;
      list.append(el('li', null, el('span', { class: ok ? 'pass' : 'fail' }, ok ? '✔ ' : '✘ '), t.name, detail ? el('span', { class: 'note' }, ' — ' + detail) : null));
    });
    main.querySelector('.statement').append(el('div', null, el('b', null, `${pass} / ${QM.tests.length} passed`)));
  }

  function labPage(topic, problemId) {
    const all = topicProblems(topic.id), sandbox = all.find((p) => p.sandbox), probs = all.filter((p) => !p.sandbox);
    const p = all.find((q) => q.id === problemId) || sandbox || probs[0];
    if (!p) { main.append(el('div', { class: 'prob-head' }, 'No presets for this lab.')); return; }
    const engine = QM.engines[p.engine];

    const sel = el('select', { class: 'param-select preset-select' },
      sandbox ? el('option', { value: sandbox.id }, '◇ ' + sandbox.title) : null,
      el('optgroup', { label: 'Problems' }, probs.map((q) => el('option', { value: q.id }, q.title + (q.sheet ? '  —  ' + q.sheet : '')))));
    sel.value = p.id;
    sel.addEventListener('change', () => { location.hash = '#/lab/' + topic.id + '/' + sel.value; });
    const host = el('div');
    const resetBtn = QM.ui.button('⟲ Reset preset', () => mountEngine());

    const head = el('div', { class: 'prob-head' },
      el('div', { class: 'lab-title' }, el('h2', null, topic.title), el('span', { class: 'chip engine' }, 'lab engine: ' + (engine ? engine.title : p.engine))),
      el('div', { class: 'note', html: topic.desc }),
      el('div', { class: 'preset-row' }, el('label', { class: 'param-label' }, 'Preset'), sel, resetBtn));
    const card = el('details', { class: 'problem-card', open: true },
      el('summary', null, p.sandbox ? p.title : 'Problem: ' + p.title),
      el('div', { class: 'prob-meta' },
        p.source ? el('a', { class: 'chip', href: '../Problems/' + encodeURIComponent(p.source), target: '_blank', title: 'Open the original PDF' }, '📄 ' + p.source) : null,
        p.sheet ? el('span', { class: 'chip' }, p.sheet) : null),
      el('div', { class: 'statement', html: p.statement || '' }));
    if (p.solution) {
      const sol = el('div', { class: 'solution', html: p.solution });
      const solBtn = QM.ui.button('Show solution', () => { sol.classList.toggle('open'); solBtn.textContent = sol.classList.contains('open') ? 'Hide solution' : 'Show solution'; });
      card.append(el('div', { class: 'prob-actions' }, solBtn), sol);
    }
    head.append(card);
    main.append(head, host);
    tex(head);
    function mountEngine() {
      unmount(); host.innerHTML = '';
      if (!engine) { host.append(el('div', { class: 'panel panel-body' }, 'Engine not found: ' + p.engine)); return; }
      try { mounted = engine.mount(host, clone(p.preset || {}), p) || null; }
      catch (e) { console.error(e); host.append(el('div', { class: 'panel panel-body' }, 'Error: ' + e.message)); }
      tex(host);
    }
    mountEngine();
  }

  /** deep copy that keeps functions (presets may contain builder callbacks) */
  function clone(x) {
    if (Array.isArray(x)) return x.map(clone);
    if (x && typeof x === 'object' && Object.getPrototypeOf(x) === Object.prototype) { const o = {}; for (const k in x) o[k] = clone(x[k]); return o; }
    return x;
  }

  function route() {
    const h = location.hash.replace(/^#/, '') || '/';
    if (h.startsWith('/p/')) { const p = R.problems.find((q) => q.id === h.slice(3)); if (p) { location.replace(labHref(p)); return; } }
    unmount();
    main.innerHTML = '';
    main.append(el('button', { class: 'btn menu-btn', onclick: () => document.getElementById('sidebar').classList.toggle('open') }, '☰ Labs'));
    const m = h.match(/^\/lab\/([^/]+)(?:\/(.+))?$/);
    if (m) { const t = R.topics.find((x) => x.id === m[1]); if (t) labPage(t, m[2]); else main.append(el('div', { class: 'prob-head' }, 'Unknown lab.')); }
    else if (h === '/sources') sourcesPage();
    else if (h === '/selftest') selftest();
    else home();
    tex(main);
    markActive();
    document.getElementById('sidebar').classList.remove('open');
    window.scrollTo(0, 0);
  }

  window.addEventListener('hashchange', route);
  function start() { buildNav(); route(); }
  if (document.readyState === 'complete') start(); else window.addEventListener('load', start);
})();
