import smtplib

from ntsparksearch.Common.Constants import Constants


class EmailSender(object):
    """Email model"""

    def __init__(self):
        self._receivers = None

    def _get_receivers(self) -> list:
        return self._receivers

    def _set_receivers(self, receivers: list) -> None:
        self._receivers = receivers

    def mail_sender(self, msg):
        server = smtplib.SMTP(host=Constants.MAIL_HOST)
        server.ehlo()
        server.starttls()
        server.login(Constants.MAIL_USER, Constants.MAIL_PASS)
        server.sendmail(Constants.MAIL_SENDER, self._receivers, msg)

    def send_email_download_initialize(self, list_of_genes_to_download: list) -> None:
        string_of_gene_ids = "\n".join(list_of_genes_to_download)

        if list_of_genes_to_download is not None:
            message_download = Constants.MSG_DOWNLOAD_INITIALIZE + string_of_gene_ids
            try:
                self.mail_sender(message_download)
            except Exception as error:
                print("error" + repr(error))

    def send_email_download_finished(self, list_of_genes_to_download_finished: list) -> None:
        string_of_gene_ids = "\n".join(list_of_genes_to_download_finished)

        if list_of_genes_to_download_finished is not None:
            message_download = Constants.MSG_DOWNLOAD_FINISHED + string_of_gene_ids
            try:
                self.mail_sender(message_download)
            except Exception as error:
                print("Exception caught at sending email operation" + repr(error))

    receivers = property(fget=_get_receivers, fset=_set_receivers)
