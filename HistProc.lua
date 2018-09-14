-- Processing histograms and remove outliers
print(" >> This script remove outliers according to histogram data << ")
-- Load input
dofile("DiffGeoStat.lua")

-- Set treshold to be consider outliers
pc = 10.0 -- limit to be consider an outlier, i. e., 1% population
cutoff = math.floor(models*pc/100.0) -- cutoff of total structures

print(" > pc = ", pc, "%")
print(" > Histogram threshold = " .. cutoff, "conformations")

-- create output file ssv
print ("@ > Writing DiffGeoSpec.ssv")
out_file = "DiffGeoSpec.ssv"
local f = assert(io.open(out_file, "w"))

-- create table residue_numbers

local residue_numbers = {}
	-- fill table
for k in pairs(hist) do
  table.insert(residue_numbers, k)
end
table.sort(residue_numbers)

--
for i = 1, #residue_numbers do
  local res, data = residue_numbers[i], hist[residue_numbers[i]]

  	-- Process curvature --
  -- create keys
  local keys = {}
  -- insert keys
  for k in pairs(data.curvature) do
    table.insert(keys, k)
  end
  table.sort(keys)

  -- Check left side
  local curvature = data.curvature[keys[1]]
  local left = curvature.pos
  for j = 1, #keys do
    local k, curvature = keys[j], data.curvature[keys[j]]
    if curvature.value > cutoff then
      left = curvature.pos
      break
    end
  end

  -- Check right side
  local curvature = data.curvature[keys[#keys]]
  local right = curvature.pos
  for j = #keys, 1, -1 do
    local k, curvature = keys[j], data.curvature[keys[j]]
    if curvature.value > cutoff then
      right = curvature.pos
      break
    end
  end
  -- Curvature interval magnitude inside the cutoff values
  local curvature_len = right - left


  -- Process torsion
  local keys = {}
  for k in pairs(data.torsion) do
    table.insert(keys, k)
  end
  table.sort(keys)

  -- Check left side
  local torsion = data.torsion[keys[1]]
  local left = torsion.pos
  for j = 1, #keys do
    local k, torsion = keys[j], data.torsion[keys[j]]
    if torsion.value > cutoff then
      left = torsion.pos
      break
    end
  end

  -- Check right side
  local torsion = data.torsion[keys[#keys]]
  local right = torsion.pos
  for j = #keys, 1, -1 do
    local k, torsion = keys[j], data.torsion[keys[j]]
    if torsion.value > cutoff then
      right = torsion.pos
      break
    end
  end
  local torsion_len = right - left

  -- Process writhing
  local keys = {}
  for k in pairs(data.writhing) do
    table.insert(keys, k)
  end
  table.sort(keys)

  -- Check left side
  local writhing = data.writhing[keys[1]]
  local left = writhing.pos
  for j = 1, #keys do
    local k, writhing = keys[j], data.writhing[keys[j]]
    if writhing.value > cutoff then
      left = writhing.pos
      break
    end
  end

  -- Check right side
  local writhing = data.writhing[keys[#keys]]
  local right = writhing.pos
  for j = #keys, 1, -1 do
    local k, writhing = keys[j], data.writhing[keys[j]]
    if writhing.value > cutoff then
      right = writhing.pos
      break
    end
  end
  local writhing_len = right - left

    -- Process arc length
  local keys = {}
  for k in pairs(data.arc_length) do
    table.insert(keys, k)
  end
  table.sort(keys)

  -- Check left side
  local arc_length = data.arc_length[keys[1]]
  local left = arc_length.pos
  for j = 1, #keys do
    local k, arc_length = keys[j], data.arc_length[keys[j]]
    if arc_length.value > cutoff then
      left = arc_length.pos
      break
    end
  end

  -- Check right side
  local arc_length = data.arc_length[keys[#keys]]
  local right = arc_length.pos
  for j = #keys, 1, -1 do
    local k, arc_length = keys[j], data.arc_length[keys[j]]
    if arc_length.value > cutoff then
      right = arc_length.pos
      break
    end
  end
  local arc_length_len = right - left

  f:write(string.format("%i %6.3f %6.3f %6.3f %6.3f\n", res, curvature_len, torsion_len, writhing_len, arc_length_len))

end

print(":: DONE ::")

f:close()
