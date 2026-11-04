// Collapsible Explore Data cards are native <details> elements. Shiny stops
// updating outputs that are hidden, and doesn't notice a <details> reopening,
// so signal it ("shown") to let the outputs inside catch up.
document.addEventListener("toggle", function (e) {
  var el = e.target;
  if (el.matches && el.matches("details.collapsible-card") && el.open && window.jQuery) {
    window.jQuery(el).trigger("shown");
  }
}, true);
