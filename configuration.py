# Default unit: milimeters (use Xm to use meters instead of mm)
configurationData = dict(
    surface=dict(
        dimensions=dict(
            # X axis
            width=333,
            # Y axis
            depth=333,
            # Z axis
            height=333,
        ),
        shape=dict(
            zeroHeight=166.5,
        )
    ),
    manufacture=dict(
        grid=dict(
            # X axis in arc lenght units (max value)
            width=5,
            # Y axis in arc lenght units (max value)
            height=5,
        ),
        mold=dict(
            puzzle=dict(
                thickness=3,
                pieces=dict(
                    # X axis in arc lenght units (max value)
                    width=50,
                    # Y axis in arc lenght units (max value)
                    height=80,
                ),
                tabs=dict(
                    # By default tab is y-axis oriented
                    file="svg/union.svg",
                    width=10,
                    thickness=2
                ),
            ),
            panels=dict(
                dimensions=dict(
                    # unit: pieces
                    width=5,   
                    # unit: pieces
                    height=7,
                ),
                leads=dict(
                    height=20,
                    thickness=5,
                    screws=2,
                    screwsDiameter=1,
                    insertNutDiameter=2,
                    insertNutDepth=3,
                    flangeSize=20,
                    flangeScrewsDiameter=2                    
                ),
            ),
        ),
        steelFrame=dict(
            steel=dict(
                thickness=0.002,
                heigh=0.02,
            ),
            maxDistanceError=0.005,
            notchWidth=0.0006,
            notchLengt=0.019,
            simulateNotchAngle=False,
            divisions=2,
        ),
        fiberGlass=dict(
            firstLayerThickness=0.005,
        ),
    )
)