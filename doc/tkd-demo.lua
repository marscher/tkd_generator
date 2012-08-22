--------------------------------------------------------------------------------
--
--   Lua - Script demonstrating the tkd-modeller
--
--   Authors: Arne Naegel, Sebastian Reiter, Martin Scherer
--
--------------------------------------------------------------------------------


-- specify geometry by tkd-modeller
dom = Domain3d()
tkdGenerator = TKDDomainGenerator(dom:grid(), dom:subset_handler())


-- specify cell
a = 10.0;
w = 30.0;
h = 1.0;
d = 0.2;

-- specify cell cluster
numRows = 1
numCols = 1
numLayers = 2

-- tkdGenerator:createSimpleTKDDomain(10,30,1, numRows, numCols, numLayers);
tkdGenerator:createSCDomain(a, w, h, d, numRows, numCols, numLayers);
tkdGenerator:setSubsetInfo("COR", "LIP", MakeVec(1, 0.5, 0, 1), MakeVec(0, 0.5, 1, 1));

SaveDomain(dom, "tkd_simple_domain.ugx")
