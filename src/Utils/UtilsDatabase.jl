function ExtractSimName(simu_name)
    match_result = match(r"d1cf.*", simu_name)
    if match_result === nothing
        return ""  # Return an empty string if no match is found
    else
        return match_result.match
    end
end

function ExtractAge(filename)
    age_str = match(r"\d+(?=\D*$)", filename).match
    parse(Float64, age_str)
end

function ExtractInfo(simu_name)
    
    age = parse(Int, match(r"(\d+)kyr", simu_name).captures[1])
    compression_percentage = parse(Int, match(r"cf(\d+)", simu_name).captures[1])*10
    rms_velocity = parse(Int, match(r"rms(\d+)nograv", simu_name).captures[1])
    
    return age, compression_percentage, rms_velocity
end

function TrimSimulationPaths(paths)
    trimmed_paths = []
    for path in paths
        trimmed_path = match(r".*256", path).match
        push!(trimmed_paths, trimmed_path)
    end
    return trimmed_paths
end
