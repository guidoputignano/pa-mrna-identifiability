// Render an HTML file holding one SVG to a vector PDF of the given size (mm) and a PNG preview.
// usage: node render_svg.js in.html out.pdf out.png widthMM heightMM
const { chromium } = require('playwright');
(async () => {
  const [inp, pdf, png, wmm, hmm] = process.argv.slice(2);
  const fs = require('fs');
  const html = fs.readFileSync(inp, 'utf8');
  const b = await chromium.launch();
  const p = await b.newPage();
  await p.setContent(html, { waitUntil: 'load' });
  await p.pdf({ path: pdf, width: `${wmm}mm`, height: `${hmm}mm`, printBackground: true });
  const W = Math.round(+wmm * 10), H = Math.round(+hmm * 10);
  const p2 = await b.newPage({ viewport: { width: W, height: H }, deviceScaleFactor: 1 });
  await p2.setContent(html.replace(`width="${wmm}mm" height="${hmm}mm"`, `width="${W}" height="${H}"`), { waitUntil: 'load' });
  await p2.screenshot({ path: png, fullPage: true });
  await b.close();
})();
