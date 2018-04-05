function hilbert2int(coords::Vector{Int64})
  nD = length(coords)
  coord_chunks = unpack_coords(coords)
  nChunks = length(coord_chunks)
  mask = 2^nD - 1
  start, finish = initial_start_end(nChunks, nD)
  #Int here?
  index_chunks = zeros(Int64,nChunks)
  for j = 0:(nChunks-1)
    i = gray_decode_travel(start, finish, mask, coord_chunks[j+1])
    index_chunks[j+1] = i
    start, finish = child_start_end(start, finish, mask, i)
  end
  return pack_index(index_chunks, nD)
end

function unpack_coords(coords::Vector{Int64})
  nD = length(coords)
  #biggest = reduce(max, coords)  # the max of all coords
  #nChunks = max(1, ceil(Int64,log(2, biggest + 1))) # max # of bits
  nChunks = max(1, ceil(Int64,log(2, maximum(coords) + 1))) # max # of bits
  return transpose_bits(coords, nChunks)
end

function pack_indexSlow(chunks::Vector{Int64}, nD::Int64)
  p = 2^nD  # Turn digits mod 2**nD back into a single number:
  #return reduce( lambda n, chunk: n * p + chunk, chunks )
  return reduce((n,chunk)->n*p+chunk, chunks)
end

function pack_index(chunks::Vector{Int64}, nD::Int64)
  p = 2^nD  # Turn digits mod 2**nD back into a single number:
  #return reduce( lambda n, chunk: n * p + chunk, chunks )
  sum = 0::Int64
  @simd for i = 1:length(chunks)
    sum *= p
    sum += chunks[i]
  end
  return sum
end

function initial_start_end(nChunks::Int64, nD::Int64)
  # This orients the largest cube so that
  # its start is the origin (0 corner), and
  # the first step is along the x axis, regardless of nD and nChunks:
  return 0, 2^mod(-nChunks-1,nD)  # in Python 0 <=  a % b  < b.
end

function gray_encode( bn::Int64)
    assert(bn >= 0)
    return bn $ fld(bn,2)
end

function gray_decode(n::Int64)
  sh = 1
  while true
    div = n >> sh
    n $= div
    if div <= 1
      return n
    end
    sh = sh << 1
  end
end

function gray_decode_travel(start::Int64, finish::Int64, mask::Int64, g::Int64)
  travel_bit = start $ finish
  modulus = mask + 1          # == 2**nBits
  rg = (g $ start)*fld(modulus ,travel_bit*2)
  return gray_decode((rg | fld(rg,modulus)) & mask)
end

## gray_encode_travel -- gray_encode given start and end using bit rotation.
#    Modified Gray code.  mask is 2**nbits - 1, the highest i value, so
#        gray_encode_travel( start, end, mask, 0 )    == start
#        gray_encode_travel( start, end, mask, mask ) == end
#        with a Gray-code-like walk in between.
#    This method takes the canonical Gray code, rotates the output word bits,
#    then xors ("^" in Python) with the start value.
#
function gray_encode_travel(start::Int64, finish::Int64, mask::Int64, i::Int64)
    travel_bit = start $ finish
    modulus = mask + 1          # == 2**nBits
    # travel_bit = 2**p, the bit we want to travel.
    # Canonical Gray code travels the top bit, 2**(nBits-1).
    # So we need to rotate by ( p - (nBits-1) ) == (p + 1) mod nBits.
    # We rotate by multiplying and dividing by powers of two:
    g = gray_encode( i ) * ( travel_bit * 2 )
    return ((g | fld(g, modulus)) & mask) $ start
end

## child_start_end( parent_start, parent_end, mask, i ) -- Get start & end for child.
#    i is the parent's step number, between 0 and mask.
#    Say that parent( i ) =
#           gray_encode_travel( parent_start, parent_end, mask, i ).
#    And child_start(i) and child_end(i) are what child_start_end()
#    should return -- the corners the child should travel between
#    while the parent is in this quadrant or child cube.
#      o  child_start( 0 ) == parent( 0 )       (start in a corner)
#      o  child_end( mask ) == parent( mask )   (end in a corner)
#      o  child_end(i) - child_start(i+1) == parent(i+1) - parent(i)
#         (when parent bit flips, same bit of child flips the opposite way)
#    Those constraints still leave choices when nD (# of bits in mask) > 2.
#    Here is how we resolve them when nD == 3 (mask == 111 binary),
#    for parent_start = 000 and parent_end = 100 (canonical Gray code):
#         i   parent(i)    child_
#         0     000        000   start(0)    = parent(0)
#                          001   end(0)                   = parent(1)
#                 ^ (flip)   v
#         1     001        000   start(1)    = parent(0)
#                          010   end(1)                   = parent(3)
#                ^          v
#         2     011        000   start(2)    = parent(0)
#                          010   end(2)                   = parent(3)
#                 v          ^
#         3     010        011   start(3)    = parent(2)
#                          111   end(3)                   = parent(5)
#               ^          v
#         4     110        011   start(4)    = parent(2)
#                          111   end(4)                   = parent(5)
#                 ^          v
#         5     111        110   start(5)    = parent(4)
#                          100   end(5)                   = parent(7)
#                v          ^
#         6     101        110   start(6)    = parent(4)
#                          100   end(6)                   = parent(7)
#                 v          ^
#         7     100        101   start(7)    = parent(6)
#                          100   end(7)                   = parent(7)
#    This pattern relies on the fact that gray_encode_travel()
#    always flips the same bit on the first, third, fifth, ... and last flip.
#    The pattern works for any nD >= 1.
#
function child_start_end(parent_start::Int64, parent_end::Int64, mask::Int64, i::Int64)
  start_i = max(0,    (i - 1) & ~1)  # next lower even number, or 0
  end_i =   min(mask, (i + 1) |  1)  # next higher odd number, or mask
  child_start = gray_encode_travel(parent_start, parent_end, mask, start_i)
  child_end   = gray_encode_travel(parent_start, parent_end, mask, end_i)
  return child_start, child_end
end

## transpose_bits --
#    Given nSrcs source ints each nDests bits long,
#    return nDests ints each nSrcs bits long.
#    Like a matrix transpose where ints are rows and bits are columns.
#    Earlier srcs become higher bits in dests;
#    earlier dests come from higher bits of srcs.

@inbounds function transpose_bits(srcs::Vector{Int64}, nDests::Int64)
  #srcs = copy(srcs)  # Make a copy we can modify safely.
  nSrcs = length(srcs)
  dests = Array{Int64}(nDests)#zeros(Int64,nDests)
  # Break srcs down least-significant bit first, shifting down:
  for j = (nDests - 1):-1:0
    # Put dests together most-significant first, shifting up:
    dest = 0
    for k = 0:(nSrcs-1)
      dest = dest*2 + mod(srcs[k+1],2)
      srcs[k+1] = fld(srcs[k+1],2)
    end
    dests[j+1] = dest
  end
  return dests
end

#     Overwrites srcs! and writes data to dests,
#apparently slower!!!!
@inbounds function transpose_bits!(srcs::Vector{Int64}, nDests::Int64, dests::Vector{Int64})
  nSrcs = length(srcs)
  # Break srcs down least-significant bit first, shifting down:
  for j = (nDests - 1):-1:0
    # Put dests together most-significant first, shifting up:
    dest = 0
    for k = 0:(nSrcs-1)
      dest = dest*2 + mod(srcs[k+1],2)
      srcs[k+1] = fld(srcs[k+1],2)
    end
    dests[j+1] = dest
  end
  return dests
end
