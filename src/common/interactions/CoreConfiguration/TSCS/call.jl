"""
Defines function call for `TSCS` on a `CoreConfiguration`
"""
function (tscs::TSCS)(design::CoreConfiguration)
    config = merge_configs(design.core, design.config)
    return tscs(config)
end