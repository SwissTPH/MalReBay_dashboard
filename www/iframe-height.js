// When the dashboard is embedded in an iframe (e.g. on the RBM site), report
// the page's content height to the parent page so it can size the iframe to
// fit, instead of clipping the content or scrolling it in a nested scrollbar.
//
// Message posted to the parent, every time the height changes:
//   { type: "malrebay:height", height: <integer, CSS pixels> }
// The parent can also ask for the current height at any time by posting
//   { type: "malrebay:request-height" }
// into the iframe.
(function () {
  if (window.self === window.top) return;   // not embedded: nothing to do

  // Lets custom.css drop rules tied to the viewport height (100vh), which
  // inside an auto-sized iframe would stop the page from ever shrinking.
  document.documentElement.classList.add("in-iframe");

  var lastHeight = -1;
  var pending = false;

  function contentHeight() {
    var body = document.body;
    if (!body) return 0;
    var style = window.getComputedStyle(body);
    return Math.ceil(body.getBoundingClientRect().height +
                     parseFloat(style.marginTop) + parseFloat(style.marginBottom));
  }

  function send(force) {
    var h = contentHeight();
    if (!force && h === lastHeight) return;
    lastHeight = h;
    window.parent.postMessage({ type: "malrebay:height", height: h }, "*");
  }

  // Batch bursts of layout changes (plots rendering, tabs switching, cards
  // opening) into one message per frame
  function schedule() {
    if (pending) return;
    pending = true;
    window.requestAnimationFrame(function () { pending = false; send(false); });
  }

  function start() {
    send(true);
    if (window.ResizeObserver) {
      new ResizeObserver(schedule).observe(document.body);
    } else {
      setInterval(schedule, 500);
    }
    // Fallbacks for changes a ResizeObserver can miss (late images, fonts)
    window.addEventListener("load", schedule);
    document.addEventListener("shiny:value", schedule);
    document.addEventListener("toggle", schedule, true);
  }

  window.addEventListener("message", function (e) {
    if (e.data && e.data.type === "malrebay:request-height") send(true);
  });

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", start);
  } else {
    start();
  }
})();
