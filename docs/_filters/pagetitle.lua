-- Give pages that have no title of their own (the quartodoc API pages) a
-- browser-tab <title> taken from their first level-1 heading, preserving case.
-- Quarto otherwise derives a lowercased `pagetitle` (mb_solve.mbsolve); this
-- overrides it. Setting `pagetitle` leaves the page body alone, so no second
-- title block appears.
function Pandoc(doc)
  if doc.meta.title then
    return nil
  end
  for _, block in ipairs(doc.blocks) do
    if block.t == "Header" and block.level == 1 then
      doc.meta.pagetitle = pandoc.utils.stringify(block.content)
      return doc
    end
  end
  return nil
end
