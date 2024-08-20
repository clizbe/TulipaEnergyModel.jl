# At the end of the file, there is a reference relating schemas and files

const schemas = (
    assets = (
        # Schema for the assets-data.csv file.
        data = OrderedDict(
            :name => "VARCHAR",                              # Name of Asset (geographical?)
            :type => "VARCHAR",                              # Producer/Consumer/Storage/Conversion
            :active => "BOOLEAN",                            # Active or decomissioned
            :group => "VARCHAR",                             # Group to which the asset belongs to (missing -> no group)
            :investable => "BOOLEAN",                        # Whether able to invest
            :investment_integer => "BOOLEAN",                # Whether investment is integer or continuous
            :investment_cost => "DOUBLE",                    # kEUR/MW/year
            :investment_limit => "DOUBLE",                   # MW (Missing -> no limit)
            :capacity => "DOUBLE",                           # MW
            :initial_capacity => "DOUBLE",                   # MW
            :peak_demand => "DOUBLE",                        # MW
            :consumer_balance_sense => "VARCHAR",            # Sense of the consumer balance constraint (default ==)
            :is_seasonal => "BOOLEAN",                       # Whether seasonal storage (e.g. hydro) or not (e.g. battery)
            :storage_inflows => "DOUBLE",                    # MWh/year
            :initial_storage_capacity => "DOUBLE",           # MWh
            :initial_storage_level => "DOUBLE",              # MWh (Missing -> free initial level)
            :energy_to_power_ratio => "DOUBLE",              # Hours
            :storage_method_energy => "BOOLEAN",             # Whether storage method is energy or not (i.e., fixed_ratio)
            :investment_cost_storage_energy => "DOUBLE",     # kEUR/MWh/year
            :investment_limit_storage_energy => "DOUBLE",    # MWh (Missing -> no limit)
            :capacity_storage_energy => "DOUBLE",            # MWh
            :investment_integer_storage_energy => "BOOLEAN", # Whether investment for storage energy is integer or continuous
            :use_binary_storage_method => "VARCHAR",         # Whether to use an extra binary variable for the storage assets to avoid charging and discharging simultaneously (missing;binary;relaxed_binary)
            :max_energy_timeframe_partition => "DOUBLE",     # MWh (Missing -> no limit)
            :min_energy_timeframe_partition => "DOUBLE",     # MWh (Missing -> no limit)
        ),

        # Schema for the assets-profiles.csv file.
        profiles_reference = OrderedDict(
            :asset => "VARCHAR",               # Asset name
            :profile_type => "VARCHAR",        # Type of profile, used to determine dataframe with source profile
            :profile_name => "VARCHAR",        # Name of profile, used to determine data inside the dataframe
        ),

        # Schema for the assets-timeframe-partitions.csv file.
        timeframe_partition = OrderedDict(
            :asset => "VARCHAR",
            :specification => "VARCHAR",
            :partition => "VARCHAR",
        ),

        # Schema for the assets-rep-periods-partitions.csv file.
        rep_periods_partition = OrderedDict(
            :asset => "VARCHAR",
            :rep_period => "INTEGER",
            :specification => "VARCHAR",
            :partition => "VARCHAR",
        ),
    ),
    groups = (
        # Schema for the groups-data.csv file.
        data = OrderedDict(
            :name => "VARCHAR",                # Name of the Group
            :invest_method => "BOOLEAN",       # true -> activate group constraints; false -> no group investment constraints
            :min_investment_limit => "DOUBLE", # MW (Missing -> no limit)
            :max_investment_limit => "DOUBLE", # MW (Missing -> no limit)
        ),
    ),
    flows = (
        # Schema for the flows-data.csv file.
        data = OrderedDict(
            :carrier => "VARCHAR",                  # (Optional?) Energy carrier
            :from_asset => "VARCHAR",               # Name of Asset
            :to_asset => "VARCHAR",                 # Name of Asset
            :active => "BOOLEAN",                   # Active or decomissioned
            :is_transport => "BOOLEAN",             # Whether a transport flow
            :investable => "BOOLEAN",               # Whether able to invest
            :investment_integer => "BOOLEAN",       # Whether investment is integer or continuous
            :variable_cost => "DOUBLE",             # kEUR/MWh
            :investment_cost => "DOUBLE",           # kEUR/MW/year
            :investment_limit => "DOUBLE",          # MW
            :capacity => "DOUBLE",                  # MW
            :initial_export_capacity => "DOUBLE",   # MW
            :initial_import_capacity => "DOUBLE",   # MW
            :efficiency => "DOUBLE",                # p.u. (per unit)
        ),

        # Schema for the flows-profiles file.
        profiles_reference = OrderedDict(
            :from_asset => "VARCHAR",          # Name of Asset
            :to_asset => "VARCHAR",            # Name of Asset
            :profile_type => "VARCHAR",        # Type of profile, used to determine dataframe with source profile
            :profile_name => "VARCHAR",        # Name of profile, used to determine data inside the dataframe
        ),

        # Schema for the flows-rep-periods-partitions.csv file.
        rep_periods_partition = OrderedDict(
            :from_asset => "VARCHAR",          # Name of Asset
            :to_asset => "VARCHAR",            # Name of Asset
            :rep_period => "INTEGER",
            :specification => "VARCHAR",
            :partition => "VARCHAR",
        ),
    ),
    timeframe = (
        # Schema for the profiles-timeframe-<type>.csv file.
        profiles_data = OrderedDict(
            :profile_name => "VARCHAR",      # Profile name
            :period => "INTEGER",            # Period
            :value => "DOUBLE",              # p.u. (per unit)
        ),
    ),

    # Schema for the rep-periods-data.csv file.
    rep_periods = (
        data = OrderedDict(
            :rep_period => "INTEGER",        # Representative period number
            :num_timesteps => "INTEGER",     # Numer of timesteps
            :resolution => "DOUBLE",         # Duration of each timestep (hours)
        ),

        # Schema for the rep-periods-mapping.csv file.
        mapping = OrderedDict(
            :period => "INTEGER",            # Period number
            :rep_period => "INTEGER",        # Representative period number
            :weight => "DOUBLE",             # Hours
        ),

        # Schema for the profiles-rep-periods-<type>.csv file.
        profiles_data = OrderedDict(
            :profile_name => "VARCHAR",  # Profile name
            :rep_period => "INTEGER",    # Representative period number
            :timestep => "INTEGER",      # Timestep number
            :value => "DOUBLE",          # p.u. (per unit)
        ),
    ),
)

const schema_per_table_name = OrderedDict(
    "assets_timeframe_partitions" => schemas.assets.timeframe_partition,
    "assets_data" => schemas.assets.data,
    "assets_timeframe_profiles" => schemas.assets.profiles_reference,
    "assets_profiles" => schemas.assets.profiles_reference,
    "assets_rep_periods_partitions" => schemas.assets.rep_periods_partition,
    "flows_data" => schemas.flows.data,
    "flows_profiles" => schemas.flows.profiles_reference,
    "flows_rep_periods_partitions" => schemas.flows.rep_periods_partition,
    "groups_data" => schemas.groups.data,
    "profiles_timeframe" => schemas.timeframe.profiles_data,
    "profiles_rep_periods" => schemas.rep_periods.profiles_data,
    "rep_periods_data" => schemas.rep_periods.data,
    "rep_periods_mapping" => schemas.rep_periods.mapping,
)
