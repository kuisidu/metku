def log_to_file(function):
    def wrapper():
        import sys
        import datetime as dt
        import inspect
        import os
        original_stdout = sys.stdout
        filename = f'{dt.datetime.now().strftime("%Y_%m_%d__%H_%M_%S")}_{function.__name__}'
        os.mkdir(f'{os.getcwd()}/{filename}')
        os.chdir(f'{os.getcwd()}/{filename}')
        print(f'writing to file: {filename}')
        with open(f'{filename}.txt', 'w', encoding='UTF-8') as file:
            sys.stdout = file

            print(inspect.getsource(sys.modules['__main__']))
            print('\n\n###################################################################################################\n\n')
            function()

            sys.stdout = original_stdout
        print(f'done writing')
    return wrapper
