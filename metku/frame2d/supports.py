
class Support:

    def __init__(self, coordinate, **kwargs):

        values = dict(
            sid=0,
            x=None,
            y=None,
            z=None,
            rx=None,
            ry=None,
            rz=None
        )
        values.update(kwargs)

        self.node = None
        self.node_id = None
        self.coordinate = coordinate
        self.sid = values['sid']


