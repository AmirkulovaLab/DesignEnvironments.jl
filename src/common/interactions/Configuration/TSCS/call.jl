"""
Defines function call for `TSCS` on a `Configuration`

# Example
```
objective(config)
```
"""
function (tscs::TSCS)(config::Configuration)
    return tscs(config.pos)
end