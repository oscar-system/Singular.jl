using DataFrames
using CSV
function savew(df::DataFrame, file::String)

  open(file, "w") do io # create a file and write with header as the file does not exist
    foreach(row -> print(io, row), CSV.RowWriter(df))
  end
end

function savea(df::DataFrame, file::String)

  open(file, "a") do io # append to file and write without header
    foreach(row -> print(io, row), CSV.RowWriter(df, writeheader = false))
  end
end

function savea(df::DataFrame, file::String, p::Int64)
  if p == 2
    open("pertubedWalk2.txt","a") do io # append to file and write without header
      foreach(row -> print(io, row), CSV.RowWriter(df, writeheader = false))
    end
  elseif p == 3
    open("pertubedWalk3.txt","a") do io # append to file and write without header
      foreach(row -> print(io, row), CSV.RowWriter(df, writeheader = false))
    end
  elseif p == 4
    open("pertubedWalk4.txt","a") do io # append to file and write without header
      foreach(row -> print(io, row), CSV.RowWriter(df, writeheader = false))
    end
  elseif p == 5
    open("pertubedWalk5.txt","a") do io # append to file and write without header
      foreach(row -> print(io, row), CSV.RowWriter(df, writeheader = false))
    end
  elseif p == 6
    open("pertubedWalk6.txt","a") do io # append to file and write without header
      foreach(row -> print(io, row), CSV.RowWriter(df, writeheader = false))
    end
  elseif p == 7
    open("pertubedWalk7.txt","a") do io # append to file and write without header
      foreach(row -> print(io, row), CSV.RowWriter(df, writeheader = false))
    end
  elseif p == 8
    open("pertubedWalk8.txt","a") do io # append to file and write without header
      foreach(row -> print(io, row), CSV.RowWriter(df, writeheader = false))
    end
  elseif p == 9
    open("pertubedWalk9.txt","a") do io # append to file and write without header
      foreach(row -> print(io, row), CSV.RowWriter(df, writeheader = false))
    end
  elseif p == 10
    open("pertubedWalk10.txt","a") do io # append to file and write without header
      foreach(row -> print(io, row), CSV.RowWriter(df, writeheader = false))
    end
  end
end
