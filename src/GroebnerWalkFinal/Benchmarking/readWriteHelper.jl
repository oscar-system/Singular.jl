using DataFrames
using CSV
function savew(df::DataFrame, file::String)

  open(file * ".txt", "w") do io # create a file and write with header as the file does not exist
    foreach(row -> print(io, row), CSV.RowWriter(df))
  end
end

function savea(df::DataFrame, file::String)

  open(file * ".txt", "a") do io # append to file and write without header
    foreach(row -> print(io, row), CSV.RowWriter(df, writeheader = false))
  end
end

function savea(df::DataFrame, file::String, p::Int64)
  open(file * "$p.txt", "a") do io # append to file and write without header
    foreach(row -> print(io, row), CSV.RowWriter(df, writeheader = false))
  end
end

function savew(df::DataFrame, file::String, p::Int64)
  open(file * "$p.txt", "w") do io # append to file and write without header
    foreach(row -> print(io, row), CSV.RowWriter(df, writeheader = false))
  end
end
