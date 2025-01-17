class Chromo(object):
    def __init__(self, **data):
        self.__dict__.update(data)  # 在__dict__属性中更新data
        self.size = len(data['data'])