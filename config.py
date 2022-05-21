# Default unit: milimeters (use Xm to use meters instead of mm)
configuration = dict(
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
            # X axis
            width=3.33,
            # Y axis
            height=3.33,
        ),
        mold=dict(
            puzzle=dict(
                thickness=10,
                pieces=dict(
                    # X axis
                    width=83.25,
                    # Y axis
                    height=55.5,
                ),
                tabs=dict(
                    # By default tab is y-axis oriented
                    file="svg/union.svg",
                    width=20,
                ),
            ),
            panels=dict(
                dimensions=dict(
                    # unit: pieces
                    width=2,   
                    # unit: pieces
                    height=3,
                ),
                leads=dict(
                    height=20,
                    thickness=5,
                    #boltsPerFace=2,
                    #boltsDiameter=0.003,
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
